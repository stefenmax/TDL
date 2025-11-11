function [D,Gamma,err,gerr] = kcpd20161017(params,varargin)
% This function is used to train a tensor dictionary using K-CPD algorithm
% This code is based on  Ron Rubinstein's KSVD package

global CODE_SPARSITY CODE_ERROR codemode
global MEM_LOW MEM_NORMAL MEM_HIGH memusage
global ompfunc ompparams exactsvd

CODE_SPARSITY = 1;
CODE_ERROR = 2;

MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;


%%%%% parse input parameters %%%%%


data = params.data;

ompparams = {'checkdict','off'};

% coding mode %

if (isfield(params,'codemode'))
  switch lower(params.codemode)
    case 'sparsity'
      codemode = CODE_SPARSITY;
      thresh = params.Tdata;
    case 'error'
      codemode = CODE_ERROR;
      thresh = params.Edata;
    otherwise
      error('Invalid coding mode specified');
  end
elseif (isfield(params,'Tdata'))
  codemode = CODE_SPARSITY;
  thresh = params.Tdata;
elseif (isfield(params,'Edata'))
  codemode = CODE_ERROR;
  thresh = params.Edata;

else
  error('Data sparse-coding target not specified');
end


% max number of atoms %

if (codemode==CODE_ERROR && isfield(params,'maxatoms'))
  ompparams{end+1} = 'maxatoms';
  ompparams{end+1} = params.maxatoms;
end


% memory usage %

if (isfield(params,'memusage'))
  switch lower(params.memusage)
    case 'low'
      memusage = MEM_LOW;
    case 'normal'
      memusage = MEM_NORMAL;
    case 'high'
      memusage = MEM_HIGH;
    otherwise
      error('Invalid memory usage mode');
  end
else
  memusage = MEM_NORMAL;
end


% iteration count %

if (isfield(params,'iternum'))
  iternum = params.iternum;
else
  iternum = 10;
end


% omp function %

if (codemode == CODE_SPARSITY)
  ompfunc = @omp;
else
  ompfunc = @omp2;
end


% status messages %

printiter = 0;
printreplaced = 1;
printerr = 0;
printgerr = 0;

verbose = 't';
msgdelta = -1;

for i = 1:length(varargin)
  if (ischar(varargin{i}))
    verbose = varargin{i};
  elseif (isnumeric(varargin{i}))
    msgdelta = varargin{i};
  else
    error('Invalid call syntax');
  end
end

for i = 1:length(verbose)
  switch lower(verbose(i))
    case 'i'
      printiter = 1;
    case 'r'
      printiter = 1;
      printreplaced = 1;
    case 't'
      printiter = 1;
      printerr = 1;
      if (isfield(params,'testdata'))
        printgerr = 1;
      end
  end
end

if (msgdelta<=0 || isempty(verbose))
  msgdelta = -1; 
end

ompparams{end+1} = 'messages';
ompparams{end+1} = msgdelta;



% compute error flag %

comperr = (nargout>=3 || printerr);


% validation flag %

testgen = 0;
if (isfield(params,'testdata'))
  testdata = params.testdata;
  if (nargout>=4 || printgerr)
    testgen = 1;
  end
end


% data norms %

tensorblocksize = params.tensorblocksize;
XtX = []; XtXg = [];
if (codemode==CODE_ERROR && memusage==MEM_HIGH)
  XtX = block_squared(data,tensorblocksize);
  if (testgen)
    XtXg = block_squared(testdata,tensorblocksize);
  end
end


% mutual incoherence limit %

if (isfield(params,'muthresh'))
  muthresh = params.muthresh;
else
  muthresh = 0.99;
end
if (muthresh < 0)
  error('invalid muthresh value, must be non-negative');
end


% exact svd computation %

exactsvd = 0;
if (isfield(params,'exact') && params.exact~=0)
  exactsvd = 1;
end


% determine dictionary size %

if (isfield(params,'initdict'))
  if (any(size(params.initdict)==1) && all(iswhole(params.initdict(:))))
    dictsize = length(params.initdict);
  else
    dictsize = size(params.initdict,2);
  end
end
if (isfield(params,'dictsize'))    % this superceedes the size determined by initdict
  dictsize = params.dictsize;
end

D = params.Dinit;   % initialize dictionary

% normalize the dictionary %

D = normatoms(D);

err = zeros(1,iternum);
gerr = zeros(1,iternum);

if (codemode == CODE_SPARSITY)
  errstr = 'RMSE';
else
  errstr = 'mean atomnum';
end


%%%%%%%%%%%%%%%%%  main loop  %%%%%%%%%%%%%%%%%

% tensors to vectors
data = vect_tensor(data,tensorblocksize, 'tensor2vect');
D = vect_tensor(D,tensorblocksize, 'tensor2vect');

for iter = 1:iternum
  
  G = [];
  if (memusage >= MEM_NORMAL)
    G = D'*D;%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  
  
  %%%%%  sparse coding  %%%%%
  
  Gamma = sparsecode(data,D,XtX,G,thresh);
  ompparams1.L = params.maxatoms;
  ompparams1.eps = thresh;
%   Gamma = mexOMP(data, D, ompparams1);
  RMSEperpix = sqrt(sum(reperror2(data,D,Gamma))/numel(data))
  
  %%%%%  dictionary update  %%%%%
  
  replaced_atoms = zeros(1,dictsize);  % mark each atom replaced by optimize_atom
  
  unused_sigs = 1:size(data,2);  % tracks the signals that were used to replace "dead" atoms.
                                 % makes sure the same signal is not selected twice
                                 
  data = vect_tensor(data,tensorblocksize, 'vect2tensor');
  D = vect_tensor(D,tensorblocksize, 'vect2tensor');
    
  p = randperm(dictsize);
  tid = timerinit('updating atoms', dictsize);
  for j = 1:dictsize
%     [D(:,p(j)),gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom_CPD(data,D,p(j),Gamma,unused_sigs,replaced_atoms,tensorblocksize);

[gamma_j, data_indices] = sprow(Gamma, p(j));
if length(data_indices) >= 0
    eval(['[D(',repmat(':,',1,length(tensorblocksize)),'p(j)),gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom_CPD(data,D,p(j),Gamma,unused_sigs,replaced_atoms,tensorblocksize);'])
    Gamma(p(j),data_indices) = gamma_j;

end

    if (msgdelta>0)
      timereta(tid, j, msgdelta);
    end
  end
  if (msgdelta>0)
    printf('updating atoms: iteration %d/%d', dictsize, dictsize);
  end
  
  data = vect_tensor(data,tensorblocksize, 'tensor2vect');
  D = vect_tensor(D,tensorblocksize, 'tensor2vect');
  
  %%%%%  compute error  %%%%%
  
  if (comperr)
    err(iter) = compute_err(D,Gamma,data);
  end
  if (testgen)
    if (memusage >= MEM_NORMAL)
      G = D'*D;
    end
    GammaG = sparsecode(testdata,D,XtXg,G,thresh);
    gerr(iter) = compute_err(D,GammaG,testdata);
  end
   
  
  %%%%%  clear dictionary  %%%%%
  
  [D,cleared_atoms] = cleardict(D,Gamma,data,muthresh,unused_sigs,replaced_atoms,tensorblocksize);
%   cleared_atoms = 0;
  
  %%%%%  print info  %%%%%
  
%   RMSEperpix = sqrt(sum(reperror2(data,D,Gamma))/numel(data))
  
  info = sprintf('Iteration %d / %d complete', iter, iternum);
  if (printerr)
    info = sprintf('%s, %s = %.4g', info, errstr, err(iter));
  end
  if (printgerr)
    info = sprintf('%s, test %s = %.4g', info, errstr, gerr(iter));
  end
  if (printreplaced)
    info = sprintf('%s, replaced %d atoms', info, sum(replaced_atoms) + cleared_atoms);
  end
  
  if (printiter)
    disp(info);
    if (msgdelta>0), disp(' '); end
  end
  
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          optimize_atom_CPD           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [atom,gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom_CPD(X,D,j,Gamma,unused_sigs,replaced_atoms,tensorblocksize)

R = 1;
options = struct;
options.Algorithm = @cpd_als;

Xvect = vect_tensor(X,tensorblocksize, 'tensor2vect');
Dvect = vect_tensor(D,tensorblocksize, 'tensor2vect');

% data samples which use the atom, and the corresponding nonzero
% coefficients in Gamma
[gamma_j, data_indices] = sprow(Gamma, j);

if (length(data_indices) < 1)
  maxsignals = 5000;
  perm = randperm(length(unused_sigs));
  perm = perm(1:min(maxsignals,end));
  E = sum((Xvect(:,unused_sigs(perm)) - Dvect*Gamma(:,unused_sigs(perm))).^2);
  [d,i] = max(E);
  eval(['atomtemp = X(',repmat(':,',1,length(tensorblocksize)),'unused_sigs(perm(i)));'])
%   [P,output] = cp_als(tensor(atomtemp),1,'init','nvecs','printitn',0,'maxiters',50,'tol',1.0e-06); % P is a rank-R K-tensor with N factors
%   Patom = ktensor(1,P.U{1},P.U{2},P.U{3});
%   Patom = tensor(Patom);
%   atom = Patom.data;
  
  [U,output] = cpd(atomtemp,R,options);
  atom = outprod(U{1}, U{2}, U{3});
  atom = atom/sqrt(sum(sum(sum(atom.^2))));
  gamma_j = zeros(size(gamma_j));
  unused_sigs = unused_sigs([1:perm(i)-1,perm(i)+1:end]);
  replaced_atoms(j) = 1;
  return;
end

maxsignals = 7000;   % max samples, cp_als is very slow if it exceeds 10K and almost stops if more than 20K
perm = randperm(length(data_indices));
perm = perm(1:min(maxsignals,end));
smallGamma = Gamma(:,data_indices(perm));
Dj = Dvect(:,j);
tempvect = Xvect(:,data_indices(perm)) - Dvect*smallGamma + Dj*gamma_j(perm);
tempblock = vect_tensor(tempvect,tensorblocksize, 'vect2tensor');
% [P,output] = cp_als(tensor(tempblock),1,'init','nvecs','printitn',0,'maxiters',50,'tol',1.0e-06); % P is a rank-R K-tensor with N factors
% Patom = ktensor(1,P.U{1},P.U{2},P.U{3});
% Patom = tensor(Patom);
% atom = Patom.data;


  [U,output] = cpd(tempblock,R,options);
  atom = outprod(U{1}, U{2}, U{3});
  atomnorm = sqrt(sum(sum(sum(atom.^2))));
  atom = atom/atomnorm;

if length(size(tempblock))== 3 && length(tensorblocksize) == 3
    gamma_j = atomnorm;
elseif length(size(tempblock))== 4 && length(tensorblocksize) == 3
%     gamma_j = reshape(P.lambda*P.U{4},1, length(data_indices));
    atomvect = vect_tensor(atom,tensorblocksize, 'tensor2vect');
    Xres = Xvect(:,data_indices)- Dvect*Gamma(:,data_indices) + Dj*gamma_j; % compute the residual of all the specific atom related signals without this atom
    gamma_j = sum(Xres.*repmat(atomvect,1,length(data_indices)));
else
    error('Please check dim!')
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             sparsecode               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Gamma = sparsecode(data,D,XtX,G,thresh)

global CODE_SPARSITY codemode
global MEM_HIGH memusage
global ompfunc ompparams

if (memusage < MEM_HIGH)
  Gamma = ompfunc(D,data,G,thresh,ompparams{:});
  
else  % memusage is high
  
  if (codemode == CODE_SPARSITY)
    Gamma = ompfunc(D'*data,G,thresh,ompparams{:});
    
  else
    Gamma = ompfunc(D'*data,XtX,G,thresh,ompparams{:});
  end
  
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             compute_err              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err = compute_err(D,Gamma,data)
  
global CODE_SPARSITY codemode

if (codemode == CODE_SPARSITY)
  err = sqrt(sum(reperror2(data,D,Gamma))/numel(data));
else
  err = nnz(Gamma)/size(data,2);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           cleardict                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [D,cleared_atoms] = cleardict(D,Gamma,X,muthresh,unused_sigs,replaced_atoms,tensorblocksize)

use_thresh = 4;  % at least this number of samples must use the atom to be kept

dictsize = size(D,2);

% compute error in blocks to conserve memory
err = zeros(1,size(X,2));
blocks = [1:3000:size(X,2) size(X,2)+1];
for i = 1:length(blocks)-1
  err(blocks(i):blocks(i+1)-1) = sum((X(:,blocks(i):blocks(i+1)-1)-D*Gamma(:,blocks(i):blocks(i+1)-1)).^2);
end

cleared_atoms = 0;
usecount = sum(abs(Gamma)>1e-7, 2);

for j = 1:dictsize
  
  % compute G(:,j)
  Gj = D'*D(:,j);
  Gj(j) = 0;
  
  % replace atom
  if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh) && ~replaced_atoms(j) )
    [y,i] = max(err(unused_sigs));
    
    dtemp = X(:,unused_sigs(i)) - D*Gamma(:,unused_sigs(i));  % find the max residual of the unused signal
    dtemptensor = vect_tensor(dtemp,tensorblocksize, 'vect2tensor');    % transformed into tensors
  	[P,output] = cp_als(tensor(dtemptensor),1,'init','nvecs','printitn',0,'maxiters',50,'tol',1.0e-06); % P is a rank-R K-tensor with N factors
    Patom = ktensor(1,P.U{1},P.U{2},P.U{3});
    Patom = tensor(Patom);
    atomtensor = Patom.data;
    D(:,j) = vect_tensor(atomtensor,tensorblocksize, 'tensor2vect');
    
    unused_sigs = unused_sigs([1:i-1,i+1:end]);
    cleared_atoms = cleared_atoms+1;
  end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            misc functions            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err2 = reperror2(X,D,Gamma)

% compute in blocks to conserve memory
err2 = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  err2(blockids) = sum((X(:,blockids) - D*Gamma(:,blockids)).^2);
end

end


function Y = block_squared(X,tensorblocksize)

% compute the sum of squared pixel values in each atom

% get the number of data
datasize = size(X);
if length(datasize) == length(tensorblocksize)
    numdat = 1; 
elseif length(datasize)-length(tensorblocksize) == 1
    numdat = datasize(end);    
else
    error('Error: Please check the dim of data.')
end
Y = zeros(1,numdat);
if numdat == 1
    Y = sum(X(:).^2);
else
    for i = 1:numdat
        eval(['xtemp = X(',repmat(':,',1,length(datasize)-1),'i);'])
        Y(i) = sum(xtemp(:).^2);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         norm atom function           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dictnorm = normatoms(dict)

% normallization of each atom

dictsize = size(dict);
atomsize = dictsize(1:end-1);
numatom = dictsize(end);
dictnorm = zeros(dictsize);
for i = 1:numatom
    eval(['atomtemp = dict(',repmat(':,',1,length(atomsize)),'i);'])
    eval(['dictnorm(',repmat(':,',1,length(atomsize)),'i)= atomtemp/sqrt(sum(atomtemp(:).^2));'])    
end

end
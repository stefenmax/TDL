% 在ompdenoise3的基础上修改得到，适用于CPD分解得到三维字典。
%
% 2014-03-30
% 2014-04-26，在ompdenoise3CPD.m的基础上修改得到，适用于普通字典
% 的各个通道图像的一致原子表达

function [y,nz] = ompdenoise3_group(params,msgdelta)
%OMPDENOISE3 OMP denoising of 3-D signals.
%  OMPDENOISE3 denoises a 3-dimensional signal using OMP denoising. The
%  function syntax is identical to OMPDENOISE, but it runs significantly
%  faster on 3-D signals. OMPDENOISE3 requires somewhat more memory than
%  OMPDENOISE (approximately the size of the input signal), so if memory is
%  limited, OMPDENOISE can be used instead.
%
%  See also OMPDENOISE.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


% parse input arguments %

x = params.x;
D = params.dict;
% blocksize = params.blocksize;
tensorblocksize = params.tensorblocksize;
dctype = params.dctype;

% % blocksize %
% if (numel(blocksize)==1)
%   blocksize = ones(1,3)*blocksize;
% end


% maxval %
if (isfield(params,'maxval'))
  maxval = params.maxval;
else
  maxval = 1;
end


% gain %
if (isfield(params,'gain'))
  gain = params.gain;
else
  gain = 1.15;
end


% maxatoms %
if (isfield(params,'maxatoms'))
  maxatoms = params.maxatoms;
else
  maxatoms = floor(prod(tensorblocksize)/2);
end


% stepsize %
if (isfield(params,'stepsize'))
  stepsize = params.stepsize;
  if (numel(stepsize)==1)
    stepsize = ones(1,3)*stepsize;
  end
else
  stepsize = ones(1,3);
end
if (any(stepsize<1))
  error('Invalid step size.');
end


% noise mode %
if (isfield(params,'noisemode'))
  switch lower(params.noisemode)
    case 'psnr'
      sigma = maxval / 10^(params.psnr/20);
    case 'sigma'
      sigma = params.sigma;
    otherwise
      error('Invalid noise mode specified');
  end
elseif (isfield(params,'sigma'))
  sigma = params.sigma;
elseif (isfield(params,'psnr'))
  sigma = maxval / 10^(params.psnr/20);
else
  error('Noise strength not specified');
end


% lambda %
if (isfield(params,'lambda'))
  lambda = params.lambda;
else
  lambda = maxval/(10*sigma);
end


% msgdelta %
if (nargin <2)
  msgdelta = 5;
end
if (msgdelta<=0)
  msgdelta = -1;
end


epsilon = sqrt(prod(tensorblocksize)) * sigma * gain;   % target error for omp


MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;

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


% compute G %

G = [];
if (memusage >= MEM_NORMAL)
  G = D'*D;
end


% verify dictionary normalization %

if (isempty(G))
  atomnorms = sum(D.*D);
else
  atomnorms = diag(G);
end
if (any(abs(atomnorms-1) > 1e-2))
  error('Dictionary columns must be normalized to unit length');
end


% denoise the signal %

ompparams.nMC = 4;   % 能量通道的个数
ompparams.eps = reshape(sqrt(prod(tensorblocksize(1:end-1))) * sigma * gain,1,ompparams.nMC);
ompparams.L = maxatoms;


nz = 0;  % count non-zeros in block representations

% the denoised signal
y = zeros(size(x));

blocknum = prod(floor((size(x)-tensorblocksize)./stepsize) + 1);
processedblocks = 0;
tid = timerinit('ompdenoise', blocknum);

for k = 1:stepsize(3):size(y,3)-tensorblocksize(3)+1
  for j = 1:stepsize(2):size(y,2)-tensorblocksize(2)+1
    
    % the current batch of blocks
    blocks = im2colstep(x(:,j:j+tensorblocksize(2)-1,k:k+tensorblocksize(3)-1),tensorblocksize,stepsize);

    % remove DC
    [blocks, dc] = remove_dc_col(blocks,dctype,tensorblocksize);

    % denoise the blocks
%     if (memusage == MEM_LOW)
%       gamma = omp2(D,blocks,[],epsilon,'maxatoms',maxatoms,'checkdict','off');
%     else
%       gamma = omp2(D'*blocks,sum(blocks.*blocks),G,epsilon,'maxatoms',maxatoms,'checkdict','off');
%     end
    gamma=myOMPerr_group(D,blocks,ompparams);
    nz = nz + nnz(gamma);
%     cleanblocks = add_dc(D*gamma, dc, 'columns');
%     cleanblocks = add_dc_col(D*gamma,dc,dctype,tensorblocksize);

blockstemp = [];
for iMC = 1:ompparams.nMC
    temp = D*gamma(:,:,iMC);
    blockstemp = [blockstemp;temp];
end

    cleanblocks = add_dc_col(blockstemp,dc,dctype,tensorblocksize);

    
    cleanvol = col2imstep(cleanblocks,[size(y,1) tensorblocksize(2:3)],tensorblocksize,stepsize);
    y(:,j:j+tensorblocksize(2)-1,k:k+tensorblocksize(3)-1) = y(:,j:j+tensorblocksize(2)-1,k:k+tensorblocksize(3)-1) + cleanvol;
    
    if (msgdelta>0)
      processedblocks = processedblocks + size(blocks,2);
      timereta(tid, processedblocks, msgdelta);
    end
    
  end
end

if (msgdelta>0)
  timereta(tid, blocknum);
end


nz = nz/blocknum;  % average number of non-zeros

% average the denoised and noisy signals
cnt = countcover(size(x),tensorblocksize,stepsize);
% y = (y+lambda*x)./(cnt + lambda);
y = y./cnt; % 2014-03-30

end

%% 计算向量数据不同方式的去直流

function [y,dc] = remove_dc_col(blocks,dctype,tensorblocksize)

% % 将每个块以列向量存储改为张量块存储
% blockst = vect_tensor(blocks,tensorblocksize, 'vect2tensor');
% [yt,dc] = remove_dcCPD(blockst,dctype,tensorblocksize);
% % 转换成向量
% y = vect_tensor(yt,tensorblocksize, 'tensor2vect');
% end
if strcmpi(dctype,'total')
    [y, dc] = remove_dc(blocks,'columns');
elseif strcmpi(dctype,'channel')
    patchpixs = prod(tensorblocksize(1:end-1));
    y = [];
    dc = [];
    for i = 1:tensorblocksize(end)
        indx = (i-1)*patchpixs+1:i*patchpixs;
        [yi, dci] = remove_dc(blocks(indx,:),'columns');
        y = [y;yi];
        dc = [dc;dci];
    end
else
    error('Error: Please check the type of removing DC.')
end
end

%% 计算向量数据不同方式的恢复直流

function blocksaddc = add_dc_col(blocks,dc,dctype,tensorblocksize)

% % 将每个块以列向量存储改为张量块存储
% blockst = vect_tensor(blocks,tensorblocksize, 'vect2tensor');
% [yt,dc] = remove_dcCPD(blockst,dctype,tensorblocksize);
% % 转换成向量
% y = vect_tensor(yt,tensorblocksize, 'tensor2vect');
% end
if strcmpi(dctype,'total')
    blocksaddc = add_dc(blocks, dc, 'columns');
elseif strcmpi(dctype,'channel')
    patchpixs = prod(tensorblocksize(1:end-1));
    blocksaddc = [];
    for i = 1:tensorblocksize(end)
        indx = (i-1)*patchpixs+1:i*patchpixs;
        blocksaddci = add_dc(blocks(indx,:), dc(i,:), 'columns');
        blocksaddc = [blocksaddc;blocksaddci];
    end
else
    error('Error: Please check the type of adding DC.')
end
end


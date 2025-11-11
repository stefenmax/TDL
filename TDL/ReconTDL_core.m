function imgReconMC = ReconTDL_core(ProjMC, projTotalWeightMC, params, Niter_OSnum, fSaveName, SaveMidResultsMode)
%
% Tensor based dictionary learning for spectral CT reconstruction
%
% Input variables:
%           'ProjMC'                            spectral sinogram, 3D
%           'projTotalWeightMC'       projection weighting, 3D
%           'params'                            parameters, including initial image, regularization 
%                                                     coefficients,parameters of dictionary, etc.
%           'Niter_OSnum'                   iteration number and subset number
%           'SaveMidResultsMode'     results saving types: 0-only save the final result; 
%                                                     >1- save the main loop results
%
% Output variables:
%           'imgReconMC'                  result
%
% Yanbo Zhang
% University of Massachusetts Lowell
% yanbozhang007@gmail.com
% 2015-07-10


%% parameters

reconsize = size(params.imgInitMC,1);
[~, NumofView, nMC] = size(ProjMC);

Niter = Niter_OSnum(1, :);          % vector of main loops numbers
NiterSum = sum(Niter);
OSnum = Niter_OSnum(2, :);      % vectior of subsets numbers
OSnumIter = [];                           % subsets numbers of each iteration
for i = 1:size(Niter_OSnum,2)
    OSnumIter = [OSnumIter, OSnum(i)*ones(1, Niter(i))];
end

lambdaDict = params.lambdaDict;

% sigmas = params.sigma;
cntpartMC = zeros(size(params.imgInitMC));
for iMC = 1:nMC
    cntpartMC(:,:,iMC) = params.cnt;
end

[params,verbose,msgdelta] = parasetCPD(params);
imgkMC = params.imgInitMC;    % initial image

%% compute the denominator of fidility part

paramsones = params;
paramsones.im = ones(reconsize);
projones = projdd(paramsones);
paramsones.reconsize = reconsize;
imgRaypixsumMC = zeros(reconsize, reconsize, nMC);
for iMC = 1:nMC
    % backprojection
    paramsones.proj = projTotalWeightMC(:, :, iMC).*projones;
    bprojones = bprojdd(paramsones);
    imgRaypixsumMC(:,:,iMC) = bprojones;
end
ratedict = sum(sum(sum(imgRaypixsumMC)))/sum(sum(sum(cntpartMC)));

%% main loop

% normalize each channel
ProjMCratio = imratioMC(ProjMC, params.imratio);
imgkMC = imratioMC(imgkMC, params.imratio);
params.projMC = ProjMCratio;

for kIter = 1:NiterSum
    OS = OSnumIter(kIter);          % subset number of current iteration
    
    % iteration of each subset
    for iOS = 1:OS
        
        % compute numerator of fidelity part
        indViewOS = iOS:OS:NumofView;
        gradFidMC = - OSSARTCoreMCwdd(params, imgkMC, projTotalWeightMC, indViewOS);
        
        % OS-SART update
        if iOS ~= OS
            imgkMC = imgkMC - NumofView/length(indViewOS)*gradFidMC./imgRaypixsumMC;
        end
        
        %         % nonnegtivity
        %         imgkMC(imgkMC<0) = 0;
        
    end
    
    % representation using tensor-dictionary
    imgkMC0 = imgkMC;
    params.x = imgkMC0;
    [imgDict, nz] = ompdenoise3CPD(params,msgdelta);
    imgkMC = imgkMC0 - (OS*gradFidMC + lambdaDict.*cntpartMC.*ratedict.*(imgkMC0 - imgDict))./(imgRaypixsumMC + ratedict.*lambdaDict.*cntpartMC);
        
    % save results
    if SaveMidResultsMode>= 1
        imgkMC =  imratioMC(imgkMC, 1./params.imratio);    
        save(strcat(fSaveName,strcat('imgkMC_',num2str(kIter),'-',num2str(iOS))),'imgkMC');
        imgkMC =  imratioMC(imgkMC, params.imratio);
    end
end

imgkMC =  imratioMC(imgkMC, 1./params.imratio);    
imgReconMC = imgkMC;
save(strcat(fSaveName, 'imgReconMC'),'imgReconMC');  % save results



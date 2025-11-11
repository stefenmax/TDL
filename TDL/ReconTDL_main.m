% This script is the tensor-dictionary based spectral CT reconstruction method
%
% Yanbo Zhang
% University of Massachusetts Lowell
% yanbozhang007@gmail.com
% 2015-07-01


%% addpath

clear
clc

% find the file folder of current file
p1 = mfilename('fullpath');
i=strfind (p1,'\');
p1=p1(1:i(end));
cd(p1)

cd('../')                             % one level up from the current folder
codepath = pwd;         
addpath(genpath(codepath))  % add all the sub-folders of this folder to the path


%% load dataset

% select image, 1-real 1; 2-real 2; 3-real 3; 4- simulation 1
ind_data = 4; 

switch ind_data
    case 1
        % the first real scanned mouse dataset
        load([codepath,'\Results\Mouse_2014-6-1-15-48-5\Alldata'],'ReconrealMC','pixelsize')
        im = flipud(ReconrealMC);        
        load('ProjMCreal_case1')
        for i = 1:size(ProjMCreal_case1,3)
            ProjMC(:,:,i) = fliplr(ProjMCreal_case1(:,:,i));
        end
%         imfilename = {'00127.tif', '00128.tif'}; 
%         projfilepath0 =  'D:\Yanbo\Data\Real Data\From_XuQiong\Mouse Data\';
%         ProjMC = projprocessed(projfilepath0, imfilename); ProjMC = projprocessed(projfilepath0, imfilename,6);
%         % system matrix path
%         filepathproj = [pwdpath,'ImgLenProj360_MouseReal512_0dot036\'];
        % load tensor based dictionary
%         load([pwdpath, 'DL_sheep\result\dict3D_real_params_2015-1-12-16-45-39.mat'])
        load([codepath,'\Results\dict_CPD_ratio_params_2015-9-7-12-38-18.mat'])
        params.projMC = ProjMC;
        params.Dsource2centr = 158;
        params.Dsource2detec = 255;
        params.NumofView = 360;
        params.NumofBin = 512;       % number of detector bins
        params.pixelsize = 18.41/512;
        params.binsize = 2*0.055;
        params.binshift = 0;               % detector shift, mm
        
    case 2
        % the second real scanned mouse dataset
        load([codepath,'\Results\MouseFBP_2015-2-2-16-54-14\Alldata'],'ReconrealMC','pixelsize')
        im = flipud(ReconrealMC);
        load('ProjMCreal_case2')
        for i = 1:size(ProjMCreal_case2,3)
            ProjMC(:,:,i) = fliplr(ProjMCreal_case2(:,:,i));
        end
%         imfilename = {'00127.tif', '00128.tif'}; 
%         projfilepath0 =  'D:\Yanbo\Data\Real Data\From_XuQiong\for Yanbo\';
%         ProjMC = projprocessed(projfilepath0, imfilename);
%         % system matrix path
%         filepathproj = [pwdpath,'ImgLenProj360_MouseReal512_0dot036\'];
        % load tensor based dictionary
%         load([pwdpath, 'DL_sheep\result\dict3D_real_params_2015-1-12-16-45-39.mat'])
        load([codepath,'\Results\dict_CPD_ratio_params_2015-9-7-12-38-18.mat'])
        params.projMC = ProjMC;
        params.Dsource2centr = 158;
        params.Dsource2detec = 255;
        params.NumofView = 360;
        params.NumofBin = 512;       % number of detector bins
        params.pixelsize = 18.41/512;
        params.binsize = 2*0.055;
        params.binshift = 0;               % detector shift, mm
        
    case 3
        disp('Please select other datasets.')
        return
                
    case 4
        % the simulated mouse dataset
%         load([codepath,'\DataPrepare\simulation\Alldata_FBP2'],'ProjEnergyNoise','ReconNoise', 'pixelsize')
        load ProjEnergyNoise
        load ReconNoise
        load pixelsize
        im = ReconNoise;
%         ProjMC = ProjEnergyNoise;  
        ProjMC = zeros(size(ProjEnergyNoise));
        for i = 1:size(ProjEnergyNoise,3)
            ProjMC(:,1,i) = ProjEnergyNoise(:,1,i);
            ProjMC(:,2:end,i) = fliplr(ProjEnergyNoise(:,2:end,i));
        end
%         % system matrix path
%         filepathproj = ['D:\Yanbo\code\TDL\DataPrepare\ImgLenProj640_Moby512_0dot075\'];
        % load tensor based dictionary
%         load(['codepath,'\Results\dict_CPD_params_2015-7-21-12-37-33.mat'])
%         load([codepath,'\Results\dict_CPD_ratio_params_2015-7-21-22-14-44.mat'])
        load('dict_CPD_ratio_params_2015-7-21-22-14-44.mat')
%         load([codepath,'\Results\dict_CPD_ratio_params_2015-8-31-11-40-2.mat'])
%         load([codepath,'\Results\dict_CPD_ratio_params_2015-8-31-17-8-51.mat'])
%         load([codepath,'\Results\dict_CPD_ratio_params_2015-9-2-12-33-39.mat'])
        params.projMC = ProjMC;
        params.Dsource2centr = 132;
        params.Dsource2detec = 180;
        params.NumofView = 640;
        params.NumofBin = 512;       % number of detector bins
        params.pixelsize = 0.075;
        params.binsize = 0.1;
        params.binshift = 0;               % detector shift, mm
        
    otherwise
        error('Please input the right image index !')
        
end

%% weighting

% if using statistical method based on the Possion model,1-use; 0- do not use
isPossionWeight = 0;
% Projection weighting, considering bad detector bin, Parker weight, 
% statistical reconstruction based Poisson noise model
projBadDetecWeightMC = ones(size(ProjMC)); 
projParkerWeightMC = ones(size(ProjMC)); 
projPossionWeightMC = exp(-ProjMC);

% Statistial reconstruction or not
if isPossionWeight ==0
    projTotalWeightMC = projBadDetecWeightMC.*projParkerWeightMC;  
else
    projTotalWeightMC = projBadDetecWeightMC.*projParkerWeightMC.*projPossionWeightMC; 
end
relaxFactor = 1;         % relaxition facter, typically set to 1
projTotalWeightMC = relaxFactor*projTotalWeightMC;      % relaxition factor

%% Dictionary setting

switch ind_data
    case 1
        % the first real scanned mouse dataset
        params.lambdaDict = 1.5;%1;%1.9;%2;%0.03*64;%0.03;
        % coding precision
        params.sigma =  0.001;%0.0008;%0.0007;%0.0015;%0.0005;%0.0015;
        % sparsity level of dictionary learning
        params.maxatoms = 8;%10;%8;6
    case 2
        % the second real scanned mouse dataset
        params.lambdaDict = 1.5;%1.9;%2;%0.03*64;%0.03;
        % coding precision
        params.sigma =  0.001;%0.0008;%0.0005;%0.0015;
        % sparsity level of dictionary learning
        params.maxatoms = 8;%10;%8;
        
    case 3
        % the simulated dataset
        params.lambdaDict = 3.2;%4.5;%2;%0.05*64;%0.05;%0.05*64;%0.03;
        % coding precision
        params.sigma = 0.002/40;%0.0022;% 0.00185;%0.0020;%0.0018;%0.0005;%0.0015;
        % sparsity level of dictionary learning
        params.maxatoms = 6;%10;%10;%8;6;4
        
    case 4
        % the simulated mouse dataset
        params.lambdaDict = 3.2;%4.5;%2;%0.05*64;%0.05;%0.05*64;%0.03;
        % coding precision
        params.sigma = 0.002;%0.0018;%0.0022;% 0.00185;%0.0020;%0.0018;%0.0005;%0.0015;
        % sparsity level of dictionary learning
        params.maxatoms = 6;%10;%10;%8;6;4
        
    otherwise
        error('Please input the right image index !')
        
end

params.reconsize = 512;
nMC = size(ProjMC, 3);
imgInitMC = zeros(params.reconsize, params.reconsize, nMC);   % initalized image

% Niter_OSnum is a two rows matrix, which specifies the numbers of main loops 
% and corresponding subsets respectively of OS
Niter_OSnum = [50; 20];%[20; 40];% [50; 40];%[10 5 3 2; 10 5 2 1]; 
NiterSum = sum(Niter_OSnum(1, :));        % number of total main loops

% parameters of k-cpd
params.tensorblocksize = [8 8 nMC];%[6 6 nChannel];
params.printitn = 5;
% type of direct currency removement: 'total' removes the total DC of
% a block, 'channel' removes the DCs of each channel in a block
params.dctype = 'channel'; 

params.dict = dict;
params.imgInitMC = imgInitMC;  
params.stepsize = [1 1 1];
params.dictsize = 1024;%1024; 512;1536;2048

params.iternum = 20;
params.memusage = 'high';

% region that need to be updated by using dictionary
dictMask = zeros(size(imgInitMC)); 

% count cover of dictionary patches overlap
cnt = countcover([size(imgInitMC, 1), size(imgInitMC, 2)],params.tensorblocksize(1:2),params.stepsize(1:2));
params.cnt = cnt;

%% Setting for results saving

% saving directory of results
filesavepath0 = [codepath, '\Results\'];
t = clock;        % current time
tstr = strcat('_',num2str(t(1)),'-',num2str(t(2)),'-',num2str(t(3)),'-',num2str(t(4)),'-',num2str(t(5)),'-',num2str(floor(t(6))));
fSaveName = strcat('Case',num2str(ind_data), tstr);     % file name
fSaveName = strcat(filesavepath0, fSaveName,'\');        
mkdir(fSaveName);       % create folder

% results saving types: 0-only save the final result; 1- save the main loop results; 2-save all results
SaveMidResultsMode = 1; 

%% Dictionary based iterative reconstruction

% params.imratio = ones(nMC,1);    % only when without normalization
imgReconMC = ReconTDL_core(ProjMC, projTotalWeightMC, params, Niter_OSnum, fSaveName, SaveMidResultsMode);

% save all the variables
save([fSaveName,'Alldata'])


% 统计迭代时，只考虑源发光子数，计算各个通道的噪声比例，来估计每个通道字典去噪的sigma
%
% 2014-03-08

function gammalog = estSigmaSIR
load('D:\Program Files\MATLAB71\Tensor and Dictionary based Spectral CT recon\code\DataPrepare\result\MouseSpectralFBPrecon_2014-2-10-2-32-14\Alldata.mat','photonNumBlank')

nMC = 4;
gamma = squeeze(photonNumBlank(1,1,:));
gamma = gamma/sum(gamma);
gammalog = log(1./gamma);    

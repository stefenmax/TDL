% 计算各个通道的噪声比例，来估计每个通道字典去噪的sigma

function gammalog = estSigma
load('D:\Program Files\MATLAB71\Tensor and Dictionary based Spectral CT recon\code\DataPrepare\result\MouseSpectralFBPrecon_2014-2-10-2-32-14\Alldata.mat','ProjEnergyPhoton')

nMC = 4;
gamma = zeros(4,1);
ProjPhotonSum = sum(ProjEnergyPhoton,3);
for iMC = 1:nMC
    photoni = ProjEnergyPhoton(:,:,iMC);
    gamma(iMC) = median(median(photoni./ProjPhotonSum));
end
gammalog = log(1./gamma);    

% 计算各个通道的噪声比例，来估计每个通道字典去噪的sigma

function gammalog = estSigma_R356(ProjEnergyPhoton)

nMC = 4;
gamma = zeros(4,1);
ProjPhotonSum = sum(ProjEnergyPhoton,3);
for iMC = 1:nMC
    photoni = ProjEnergyPhoton(:,:,iMC);
    gamma(iMC) = median(median(photoni./ProjPhotonSum));
end
gammalog = log(1./gamma);    

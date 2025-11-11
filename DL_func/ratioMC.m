% 计算各个通道投影数据能量比例，以便使各个通道图像归一化（相同能量）
%
% 2014-07-07

function rMC = ratioMC(ProjMC)

projsize = size(ProjMC);
nMC = projsize(end);

proj2sumMC = zeros(1,nMC);  % 各个通道投影总能量
rMC = zeros(1,nMC); 
for iMC = 1:nMC
    proj2sumMC(iMC) = sum(sum(ProjMC(:,:,iMC).^2));       
end

proj2sumtotal = sum(proj2sumMC);    % 所有投影数据总能量

for iMC = 1:nMC
    rMC(iMC) = sqrt(1/(nMC*proj2sumMC(iMC)/proj2sumtotal));
end

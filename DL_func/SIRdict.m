% 统计迭代重建，考虑了投影权重和不同的子集个数及迭代次数
%
% 2013-12-09_9:10
% 2013-12-13，在SIR的基础上加入字典重建

function [imgRecon, normFidDict] = SIRdict(Proj, projTotalWeight, params, Niter_OSnum, pixelsize,  filepathproj, fSaveName, SaveMidResultsMode)

% params：包含初始图像，mask, 字典信息，以及字典约束项的系数等等

%% 参数计算

reconsize = size(params.imgInit,1);
[NumofRou, NumofView] = size(Proj);

Niter = Niter_OSnum(1, :);          % 大迭代次数向量
NiterSum = sum(Niter);
OSnum = Niter_OSnum(2, :);      % 子集个数向量
OSnumIter = [];                              % 每次迭代的子集个数
for i = 1:size(Niter_OSnum,2)
    OSnumIter = [OSnumIter, OSnum(i)*ones(1, Niter(i))];    
end

beta = params.beta;                     % 字典约束项的系数
imgInit = params.imgInit;             % 初始图像
dictMask = params.mask;          % 需要进行字典约束的ROI图像区域
ROIx = params.ROIx;
ROIy = params.ROIy;
cnt = params.cnt;                         % 整个图像的countcover
cntpart = cnt(ROIx, ROIy);            % 图像ROI区域的countcover
[params,verbose,msgdelta] = paraset(params); % 计算其他参数（这里用处不大）

f_DATA = reshape(imgInit',1,reconsize^2);    % 初始值选择为先验图像
imgROItemp = imgInit(ROIx, ROIy);
normFidDict = [];

%% 计算投影保真项贡献的分母

% 计算投影保真项贡献的分母
raypixsum = zeros(1,reconsize^2);                       % 计算每个像素点迭代的分母
for iView = 1:NumofView
    load(strcat(filepathproj,'ImgLen_',num2str(iView)));        % 加载该切片数据
    eval(strcat('ImgLen = ImgLen_',num2str(iView),';'));    % 将加载的该切片数据赋给统一变量ImgLen
    ImgLen(ImgLen(:,:,1)==0 ) = 1;                % 将所有未穿过的像素（为0）编号设为1，由于对应长度为0，故不影响结果
    for j = 1:NumofRou                                  % 对射线穿过的像素进行迭代
        raypixsum(ImgLen(j,:,1)) = raypixsum(ImgLen(j,:,1)) + projTotalWeight(j,iView)*sum(ImgLen(j,:,2))*ImgLen(j,:,2);
    end
    eval(strcat('clear ImgLen_',num2str(iView),';'));       % 将加载的该切片数据清除
end
raypixsum(1) = raypixsum(2);        % 由于ImgLen进行计算时对第一个像素估计不准，所以用其相邻像素代替，这种进行对结果无影响    
imgRaypixsum = reshape(raypixsum, reconsize, reconsize);
imgRaypixsum = imgRaypixsum';

% 迭代重建的分母
grad2 = imgRaypixsum + 2 * beta * cnt; 

%% 主循环，大迭代

for kIter = 1:NiterSum
    OS = OSnumIter(kIter);          % 本次迭代的子集个数
    
    % 每个子集迭代
    for iOS = 1:OS
        
        % 计算投影保真项部分的分子
        pixtemp = zeros(1,reconsize^2);          	% 计算每个像素点迭代的分子
        for iView = iOS:OS:NumofView            
            load(strcat(filepathproj,'ImgLen_',num2str(iView)));        % 加载该切片数据
            eval(strcat('ImgLen = ImgLen_',num2str(iView),';'));    % 将加载的该切片数据赋给统一变量ImgLen
            ImgLen(ImgLen(:,:,1)==0 ) = 1;                % 将所有未穿过的像素（为0）编号设为1，由于对应长度为0，故不影响结果
            
            for j = 1:NumofRou                    	% 对射线穿过的像素进行迭代                                
                pdifftemp = -(Proj(j,iView) - sum(ImgLen(j,:,2).*f_DATA(ImgLen(j,:,1)))*pixelsize/10);
                pixtemp(ImgLen(j,:,1)) = pixtemp(ImgLen(j,:,1)) + projTotalWeight(j,iView)*ImgLen(j,:,2)*pdifftemp/pixelsize*10;                
            end
            
            eval(strcat('clear ImgLen_',num2str(iView),';'));       % 将加载的该切片数据清除            
        end
        gradFid = reshape(pixtemp, reconsize, reconsize);
        gradFid = gradFid';
        
        % 计算字典约束项部分的分子
        params.x = imgROItemp;
        [imgDict, nz] = ompdenoise2(params,msgdelta);
        gradDict = zeros(size(imgInit));
        gradDict(ROIx, ROIy) = (imgROItemp - imgDict).*cntpart;               
        
        % 更新数据
        imgUpdate = (reshape(f_DATA, reconsize, reconsize))' - (OS*gradFid + 2*beta* gradDict)./grad2;
        f_DATA = reshape(imgUpdate',1,reconsize^2);
        imgROItemp = imgUpdate(ROIx, ROIy);
        
        % 计算本次迭代保真项和字典项更新强度
        normFid = norm(OS*gradFid ./ grad2, 'fro');
        normDict = norm( 2*beta* gradDict./grad2, 'fro');
        normFidDict = [normFidDict; normFid, normDict];
        
        % 保存结果
        if SaveMidResultsMode == 2
            save(strcat(fSaveName,strcat('f_DATA_',num2str(kIter),'-',num2str(iOS))),'f_DATA');  % 保存全部中间结果
        elseif SaveMidResultsMode == 1 && iOS == OS
            save(strcat(fSaveName,strcat('f_DATA_',num2str(kIter),'-',num2str(iOS))),'f_DATA');  % 只保存每次大迭代后的中间结果
        end
    end
    
end

imgRecon = reshape(f_DATA, reconsize, reconsize);
imgRecon = imgRecon';
save(strcat(fSaveName, 'imgRecon'),'imgRecon');  % 保存结果

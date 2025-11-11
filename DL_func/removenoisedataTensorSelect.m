% 从体数据块中取出各个通道方差之和最大的几个数据块
% 
% 2014-07-06

function params = removenoisedataTensorSelect(params, Num)

% datatemp = params.data;
% % datatemp2 = [];
% % bwdata = zeros(1,size(datatemp,1));
% % for i = 1:size(datatemp,2)
% %     tmp = datatemp(:,i);
% %     if std(tmp)>stdVal;
% %         datatemp2 = [datatemp2 tmp];
% %     end
% % end
% % params.data = datatemp2;
% 
% 
% stddata = std(datatemp);
% ind = stddata>stdVal;
% params.data = datatemp(:,ind);

%% 考虑各个通道的方差

arrysize = size(params.data);
dataNum = arrysize(end);  % 数据个数
if dataNum <= Num
    return
end

% 各通道去均值后重排成矩阵
tensorblocksize = params.tensorblocksize;
if strcmpi(params.dctype,'channel')
    datatemp = params.data;
else
    datatemp = remove_dcCPD(params.data,'channel',tensorblocksize);
end
datatemp2 = reshape(datatemp,prod(tensorblocksize),dataNum);

% 计算方差并排序
stddata = std(datatemp2);
stddatasort = sort(stddata,'descend');
stdVal = stddatasort(Num);
ind = stddata>=stdVal;
params.data = params.data(:,:,:,ind);






% 该函数实现对向量数据和高纬数据之间的转换
%
% 2014-03-26

function outdata = vect_tensor(data,tensorblocksize, changetype)

% input:
%   data: 输入的数据阵列（可能是高纬数据或向量）
%   tensorblocksize: 每个张量小块的size大小
%   changetype: char类型，为'vect2tensor'或'tensor2vect'
% output:
%   data: 输出的数据阵列（可能是向量或高纬数据，与输入数据形式相反）

numdata = numel(data);              % 元素总个数
numeletensor = prod(tensorblocksize);    % 每个张量的小块的元素个数
numtensor = numdata/numeletensor;   % 张量的小块的个数

if mod(numdata, numeletensor) ~= 0
    error('Please check the size of tensor!')
end

if strcmpi(changetype,'tensor2vect')
    outdata = reshape(data, [numeletensor,numtensor]);    
elseif strcmpi(changetype,'vect2tensor')
    outdata = reshape(data, [tensorblocksize,numtensor]);    
else
    error('Please input right changetype!')
end
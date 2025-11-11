% 在myOMPerr.m的基础上修改得到，用于对多通道的数据块进行联合的OMP
% 操作，即各个共同选择相同的若干原子进行表达。
%
% 2014-04-25
% 20

function [A]=myOMPerr_group_sep(D,X,ompparams); 

% 注意：这里D为一般的K-SVD得到的字典，X为转化为一维向量的
% 各个多通道数据块，所以这里每个信号的长度等于每个原子的长度
% 乘以通道个数。ompparams中增加了一个参数nMC，为通道个数，
% 而原有的参数eps则由原来的一个数值变成了nMC个数值，表示各
% 个通道相应的像素残差。返回的系数矩阵A也由原来的二维变为三
% 维，增加的第三维是表示通道数。

%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  errorGoal - the maximal allowed representation error for
%                  each siganl.
% output arguments: A - sparse coefficient matrix.
%=============================================


[ntotal,P]=size(X); % 各个通道信号大小之和ntotal，个数P
[n,K]=size(D); %一个通道信号大小n，字典个数K
% E2 = errorGoal^2*n; % 每个信号的目标误差能量
% maxNumCoef = 5;%n/2; % 表达一个信号最多的非零系数个数不能超过原信号大小的一半
nMC = ompparams.nMC;   % 能量通道的个数
E2 = ompparams.eps;
maxNumCoef = ompparams.L;
% 维度判断
if ntotal/n ~= nMC || length(E2) ~= nMC
    error('Please check dims of the argues!!!')
end

A = zeros(size(D,2),size(X,2), nMC); % 创建一个稀疏矩阵，初始为全0
info = 0;
for k=1:1:P, %按信号顺序
    a=[];
    for iMC = 1:nMC
        x=X((iMC-1)*n+1:iMC*n,k); %第k个信号的第iMC个通道
        
        residual=x;
        indx = [];
        a = [];
        currResNorm2 = sum(residual.^2, 1); % 计算各个通道
        j = 0;
        while sqrt(currResNorm2)-E2(iMC)>0 & j < maxNumCoef, % 对第k个信号只要某个通道残差大于相应阈值或编码最多更新maxNumCoef次
            j = j+1;
            proj=D'*residual;
            pos=find(abs(proj)==max(abs(proj))); % 找到proj中绝对值最大的一项的索引
            pos=pos(1);
            indx(j)=pos; %
            a=pinv(D(:,indx(1:j)))*x; % A=pinv(B),A*B*A=A,B*A*B=B, A*B is Hermitian。
            residual=x-D(:,indx(1:j))*a;
            currResNorm2 = sum(residual.^2, 1);
        end;
        if (length(indx)>0)
            A(indx,k,iMC)=a;
        end
    end
   
%    if floor(k/100)>info
%        info = k/100
%    end
end;
return;
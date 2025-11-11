function [A]=myOMPerr(D,X,ompparams); 
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  errorGoal - the maximal allowed representation error for
%                  each siganl.
% output arguments: A - sparse coefficient matrix.
%=============================================


[n,P]=size(X); % 信号大小n，个数P
[n,K]=size(D); % 字典个数K
% E2 = errorGoal^2*n; % 每个信号的目标误差能量
% maxNumCoef = 5;%n/2; % 表达一个信号最多的非零系数个数不能超过原信号大小的一半
E2 = ompparams.eps;
maxNumCoef = ompparams.L;
A = sparse(size(D,2),size(X,2)); % 创建一个稀疏矩阵，初始为全0
info = 0;
for k=1:1:P, %按信号顺序
    a=[];
    x=X(:,k); %第k个信号
    residual=x;
	indx = [];
	a = [];
	currResNorm2 = sum(residual.^2);
	j = 0;
    while currResNorm2>E2 & j < maxNumCoef, % 对第k个信号的编码最多更新maxNumCoef次
		j = j+1;
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj))); % 找到proj中绝对值最大的一项的索引
        pos=pos(1);
        indx(j)=pos; % 
        a=pinv(D(:,indx(1:j)))*x; % A=pinv(B),A*B*A=A,B*A*B=B, A*B is Hermitian。
        residual=x-D(:,indx(1:j))*a;
		currResNorm2 = sum(residual.^2);
   end;
   if (length(indx)>0)
       A(indx,k)=a;
   end
%    if floor(k/100)>info
%        info = k/100
%    end
end;
return;
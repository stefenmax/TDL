function [y,dc] = remove_dcCPD(x,dctype,tensorblocksize)
%REMOVE_DC Remove DC channel from signals.
%   [Y,DC] = REMOVE_DC(X) removes the DC channel (i.e. the mean) from the
%   specified (possibly multi-dimensional) signal X. Y is the DC-free
%   signal and is the same size as X. DC is a scalar containing the mean of
%   the signal X.
%
%   [Y,DC] = REMOVE_DC(X,'columns') where X is a 2D matrix, treats the
%   columns of X as a set of 1D signals, removing the DC channel from each
%   one individually. Y is the same size as X and contains the DC-free
%   signals. DC is a row vector of length size(X,2) containing the means of
%   the signals in X.
%
%   See also ADD_DC.


% x = params.data;
% dctype = params.dctype;
% tensorblocksize = params.tensorblocksize;

% 判断数据的个数
datasize = size(x);
if length(datasize) == length(tensorblocksize)
    numdat = 1;                 % 数据的个数
    nChannel = datasize(end);   % 每个数据的通道数
elseif length(datasize)-length(tensorblocksize) == 1
    numdat = datasize(end);    
    nChannel = datasize(end-1);  % 每个数据的通道数
else
    error('Error: Please check the dim of data.')
end

% 两种去直流方式：1：'total',整体直流为零；2：'channel',各个通道的数据均为零
if strcmpi(dctype,'total')
    if numdat == 1
        dc = mean(x(:));
        y = x - dc;
    else
        dc = zeros(1,numdat);
        y = zeros(datasize);
        for i = 1:numdat
            eval(['xtemp = x(',repmat(':,',1,length(datasize)-1),'i);'])
            dctemp = mean(xtemp(:));
            dc(i) = dctemp;
            eval(['y(',repmat(':,',1,length(datasize)-1),'i) = xtemp - dctemp;'])
        end
    end
    
elseif strcmpi(dctype,'channel')
    if numdat == 1
        dc = zeros(nChannel,1);
        for iMC = 1:nChannel
            eval(['xtemp = x(',repmat(':,',1,length(datasize)-1),'iMC);'])
            dc(iMC) = mean(xtemp(:));
            eval(['y(',repmat(':,',1,length(datasize)-1),'iMC) = xtemp - dc(iMC);'])
        end        
    else
        dc = zeros(nChannel,numdat);
        y = zeros(datasize);
        for i = 1:numdat
            eval(['xtemp0 = x(',repmat(':,',1,length(datasize)-1),'i);'])
            for iMC = 1:nChannel
                eval(['xtemp = xtemp0(',repmat(':,',1,length(datasize)-2),'iMC);'])
                dc(iMC,i) = mean(xtemp(:));
                eval(['y(',repmat(':,',1,length(datasize)-2),'iMC,i) = xtemp - dc(iMC,i);'])
            end
        end    
    
    end
    
else
    error('Error: Please check the type of removing DC.')
end

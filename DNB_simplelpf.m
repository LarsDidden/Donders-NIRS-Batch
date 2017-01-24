function [y]=DNB_simplelpf(Y)
% Simple filtering of the data using Matlabs 'butter' function.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

    index_fc = find(INFO.filt.LPF == ',');
    if isempty(index_fc) == 1
        fc = 0.5;
    else 
        fc = str2num(INFO.filt.LPF(index_fc+1:end));
    end
    [B,A] = butter(4,fc/(INFO.conv.fs/2),'low');
    for i=1:INFO.conv.total_ch
        y(:,i)=filtfilt(B,A,Y(:,i));
    end
    

end


function [y]=DNB_simplepp(Y)
% Simple detrending of the data using Matlabs 'detrend' function.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO

    for i=1:size(Y,2)
        y(:,i)=detrend(Y(:,i));
    end
    
end
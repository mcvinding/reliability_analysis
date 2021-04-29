function [alpha] = kripAlphaN2fast(dat)
% Calculate Krippendorff's alpha N=2 fast method
% Use as:
%   alpha = kripAlpha(X)
% Where 
%   dat:    N observers x M observations. For time series M = t

tic
if size(dat, 1) > 2
    error('Error: the \"N2fast\" method only works for N=2 observers. This dataset has %i.', size(dat, 1))
else
    fprintf('This dataset has %i observers and %i observations.\n',size(dat, 1) , size(dat, 2) )
end

if any(isnan(dat(:)))
    error('Error: the \"N2fast\" method is invalid for data with missing cases.\nThis data has %i missing cases.', sum(isnan(dat(:))))
end

% Get variables
allvals = unique(dat(:));
n__ = length(dat(:));

% Nominator
Zu = diff(dat).^2;                  %Interval

% denominator
dE = hist(dat(:), allvals)';        % Expected count

Zncnk = zeros(1, length(allvals));
for ii = 1:length(allvals)-1
    c = allvals(ii);                % Real value
    kvals = allvals(ii+1:end);
    deltas = (kvals-c).^2;
        
    n_c = dE(ii);
    n_k = dE(ii+1:end);
    Zncnk(ii) = sum(n_c*n_k.*deltas);
end


Do = sum(Zu);
De = sum(Zncnk) / (n__-1);

alpha = 1 - (Do/De);

dt = toc;
fprintf('Calculation done (%.3f sec).\n', dt)
%END
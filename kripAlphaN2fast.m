function [alpha, cfg] = kripAlphaN2fast(dat, method)
% Calculate Krippendorff's alpha N=2 fast method
% Use as:
%   [alpha, cfg] = kripAlpha(dat, method)
% Where 
%   dat:    N observers x M observations. For time series M = t. Note that
%           this method only works for N=2
%   method:     The method for calculating the error (i.e. delta) for
%               Krippendorpf's Alpha. Options can be:
%                 'n2fast_interval': faster computation of alpha for INTERVAL
%                   data with N=2 observers.
%                 'n2fast_nominal': faster computation of alpha for NOMINAL 
%                   data with N=2 observers.

if nargin < 2
    method = 'n2fast';
end
if size(dat, 1) > 2
    error('Error: the \"N2fast\" method only works for N=2 observers. This dataset has %i.', size(dat, 1))
else
    fprintf('This dataset has %i observers and %i observations.\n',size(dat, 1) , size(dat, 2) )
end

% if any(isnan(dat(:)))
%     error('Error: the \"N2fast\" method is invalid for data with missing cases.\nThis data has %i missing cases.', sum(isnan(dat(:))))
% end

% Get variables
allvals = unique(dat(~isnan(dat)));
n__ = length(dat(~isnan(dat)));

% Nominator
switch method
    case {'n2fast_interval', 'n2fast'}
        Zd = diff(dat).^2;                  % Interval
    case 'n2fast_nominal'
        Zd = dat(1,:)==dat(2,:);            % Nominal
end
Zu = sum(Zd(~isnan(Zd)));

% denominator
dE = hist(dat(:), allvals)';        % Expected count

% Zncnk = zeros(1, length(allvals));
Zncnk = 0;
parfor ii = 1:length(allvals)-1
    c = allvals(ii);                % Real value
    kvals = allvals(ii+1:end);
    switch method
        case {'n2fast_interval', 'n2fast'}
            deltas = (kvals-c).^2;
        case 'n2fast_nominal'
            deltas = ~(kvals==c);
    end
    n_c = dE(ii);
    n_k = dE(ii+1:end);
    Zncnk = Zncnk + sum(n_c*n_k.*deltas);
end

alpha = 1 - (Zu/Zncnk) * (n__-1);

% Variables for bootstrapping
cfg.n__     = n__;
cfg.mu      = sum(~isnan(dat));
cfg.dE      = dE;
cfg.allvals = allvals;
cfg.scale   = 'n2fast';

% dt = toc;
% fprintf('Calculation done (%.3f sec).\n', dt)
%END
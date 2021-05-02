function [alphaP] = alphaprime(dat, res, makeplot)
% Calculate Krippendorff's alpha  based on histogram density approximation 
% for faster computation on large datasets.
%
% Use as:
%   [alphaP] = alphaprime(dat, res, intpmtd, makeplot))
% Where 
%   dat:      N observers x M observations. For time series M = t
%   res:      Resolution of bins in pct of full data range (default = 0.01)
%   makeplot: Plot the data, 1 or 0 (default = 0).


%%
tic

if ~exist('res', 'var')
    res = 0.01;     % 1% bin resolution
elseif res <= 0 || res >= 0.5
    error('Resultion must be larger than 0 and below 0.5')
end

% if ~exist('method', 'var')
%     intpmtd = 'pchip';     % shape-preserving piecewise cubic interpolation
% elseif ~any(strcmp(intpmtd, {'linear','nearest','next','previous','spline','pchip','cubic','v5cubic'}))
%     error('Unknown method. See documentation for INTERP1 for valid options');
% end

if ~exist('makeplot', 'var')
    makeplot = 0;
end

fprintf('This dataset has %i observers and %i observations.\n',size(dat, 1) , size(dat, 2) )
    
% Init 
N = round(1/res);               % Number of bins
allvals = unique(dat(:));       % All unique values

% Find histogram probability function over time
% Get joint density function over time
gridx1 = linspace(min(allvals), max(allvals), N)';   % Data axis. Evenly spaces mased on max/min of data
tdim = size(dat, 2);                             % "Time" axis (copy from real data)

% X = interp1(gridx1, gridx1, dat(:), 'spline');
Y = repmat(1:size(dat, 2), size(dat, 1), 1);
XY = [dat(:), Y(:)];
dO = hist3(XY, {gridx1, 1:tdim}) ;
dE = sum(dO, 2);

n__ = sum(sum(dO));              %length(dat(:));
nu_ = sum(dO, 1);                 %size(dat,1);

% Optional plot
if makeplot
    figure;
    subplot(1,20,1:3); plot(dE, gridx1); axis tight; hold on
    subplot(1,20,5:20); imagesc(1:size(dat,2),gridx1,dO); colorbar
    set(gca,'YDir','normal')
end

% nominator
Zu = zeros(1, tdim);
for tt = 1:tdim
    if ~(nu_(tt) > 1)
        continue
    end
    % Find values
    Znucnuk = zeros(1, length(gridx1)-1);  % Initiate with zeroes.
    vidx = find(dO(:,tt));                  % Find only index with observed values, all else will be zero anyway
    
    for ii = 1:length(vidx)-1
        idx = vidx(ii);
        c = gridx1(idx);
        kvals = gridx1(vidx(ii+1:end));
        deltas = (kvals-c).^2;

        nuc = dO(idx,tt);
        nuk = dO(vidx(ii+1:end), tt);
        Znucnuk(idx) = sum(nuc*nuk.*deltas);    
    end   
    Zu(tt) = sum(Znucnuk./(nu_(tt)-1));
end

% denominator
Zncnk = zeros(1, length(gridx1));
for ii = 1:length(gridx1)-1
    c = gridx1(ii);
    kvals = gridx1(ii+1:end);
    deltas = (kvals-c).^2;
    
    n_c = dE(ii);
    n_k = dE(ii+1:end);
    Zncnk(ii) = sum(n_c*n_k.*deltas);
end

Do = sum(Zu);
De = sum(Zncnk) / (n__-1);

alphaP = 1 - (Do/De);

dt = toc;
fprintf('Calculation done (%.3f sec).\n', dt)
%END       
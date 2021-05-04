function [alphaP] = alphaprime(dat, res)
% Calculate Krippendorff's alpha  based on histogram density approximation 
% for faster computation on large datasets.
%
% Use as:
%   [alphaP] = alphaprime(dat, res, intpmtd, makeplot))
% Where 
%   dat:      N observers x M observations. For time series M = t
%   res:      Resolution of bins in pct of full data range (default = 0.01)

%%
tic

if ~exist('res', 'var')
    res = 0.01;     % 1% bin resolution
elseif res <= 0 || res >= 0.5
    error('Resultion must be larger than 0 and below 0.5')
end

fprintf('This dataset has %i observers and %i observations.\n',size(dat, 1) , size(dat, 2) )
    
% Init 
N = round(1/res);                   % Number of bins
allvals = unique(dat(:));           % All unique values
Nk = linspace(0+res, 100-res, N);

% Find histogram probability function over time
gridx = prctile(allvals, Nk)';      % Data axis. Find percentiles based on resolution.
tdim = size(dat, 2);                % "Time" axis (copy from real data)

% Make histograms
Y = repmat(1:size(dat, 2), size(dat, 1), 1);
XY = [dat(:), Y(:)];
dO = hist3(XY, {gridx, 1:tdim});    % Oberservation
dE = sum(dO, 2);                    % Expected

n__ = sum(sum(dO));              %length(dat(:));
nu_ = sum(dO, 1);                %size(dat,1);

% nominator
Zu = zeros(1, tdim);
for tt = 1:tdim
    if ~(nu_(tt) > 1)
        continue
    end
    % Find values
    Znucnuk = zeros(1, length(gridx)-1);    % Initiate with zeroes.
    vidx = find(dO(:,tt));                  % Find only index with observed values, all else will be zero anyway
    
    for ii = 1:length(vidx)-1
        idx = vidx(ii);
        c = gridx(idx);
        kvals = gridx(vidx(ii+1:end));
        deltas = (kvals-c).^2;

        nuc = dO(idx,tt);
        nuk = dO(vidx(ii+1:end), tt);
        Znucnuk(idx) = sum(nuc*nuk.*deltas);    
    end   
    Zu(tt) = sum(Znucnuk./(nu_(tt)-1));
end

% denominator
Zncnk = zeros(1, length(gridx));
for ii = 1:length(gridx)-1
    c = gridx(ii);
    kvals = gridx(ii+1:end);
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
function [alphaP] = alphaprime(dat, res, intpmtd, makeplot)
% Calculate Krippendorff's alpha  based on interpolated density approximation 
% for faster computation on large datasets.
%
% Use as:
%   [alphaP] = alphaprime(dat, res, intpmtd, makeplot))
% Where 
%   dat:      N observers x M observations. For time series M = t
%   res:      Resolution of bins in pct of full data range (default = 0.01)
%   intpmtd:  The method for interpolating values (see INTERP1 for valid
%             options). Default = 'pchip'.
%   makeplot: Plot the data, 1 or 0 (default = 0).


%%
if ~exist('res', 'var')
    res = 0.01;     % 1% bin resolution
elseif res <= 0 || res >= 0.5
    error('Resultion must be larger than 0 and below 0.5')
end

if ~exist('method', 'var')
    intpmtd = 'pchip';     % shape-preserving piecewise cubic interpolation
elseif ~any(strcmp(intpmtd, {'linear','nearest','next','previous','spline','pchip','cubic','v5cubic'}))
    error('Unknown method. See documentation for INTERP1 for valid options');
end

if ~exist('makeplot', 'var')
    makeplot = 0;
end

fprintf('This dataset has %i observers and %i observations.\n',size(dat, 1) , size(dat, 2) )
    
% Init 
N = round(1/res);               % Number of bins
allvals = unique(dat(:));       % All unique values

% Find histogram probability function over time
% Get joint density function over time
gridx1 = linspace(min(allvals), max(allvals), N);   % Data axis. Evenly spaces mased on max/min of data
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

% Nominator
Zu = zeros(1, tdim);
for tt = 1:tdim
    Vq = interp1(gridx1, dO(:,tt), allvals, intpmtd); 
    
    % normalize
    % divide by zero? No there must be a minimum of one valid value!
    c = sum(Vq)/nu_(tt);
    Vqall = Vq/c;
    
    % find test
    ZnucnukT = zeros(1, length(allvals)-1);  % Initiate with zeroes.
    vidx = find(Vqall);                      % Find only index with observed values, all else will be zero anyway
      
    for ii = 1:length(vidx)  
        idx = vidx(ii);
        c = allvals(idx);
        kvals = allvals(idx+1:end);
        deltas = (kvals-c).^2;

        nuc = Vqall(idx);
        nuk = Vqall(idx+1:end);
        ZnucnukT(idx) = sum(nuc*nuk.*deltas);  
    end
    Zu(tt) = sum(ZnucnukT./(nu_(tt)-1));
end

% denominator
V = interp1(gridx1, dE, allvals, intpmtd); 

% normalize
Zncnk = zeros(1, length(allvals));
c = sum(V)/n__;
Vall = V/c;
  
for ii = 1:length(allvals)-1
    c = allvals(ii);            % Real value
    kvals = allvals(ii+1:end);
    deltas = (kvals-c).^2;
        
    n_c = Vall(ii);
    n_k = Vall(ii+1:end);
    Zncnk(ii) = sum(n_c*n_k.*deltas);
end


Do = sum(Zu);
De = sum(Zncnk) / (n__-1);

alphaP = 1 - (Do/De);

%END       
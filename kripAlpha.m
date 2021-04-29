function [alpha] = kripAlpha(dat, scale, makeplot)
% Calculate Krippendorff's alpha using original approach.
% Use as:
%   alpha = kripAlpha(X, method, makeplot)
% Where 
%   dat:      N observers x M observations. For time series M = t
%   scale:    The method for calculating the error (i.e. delta) for
%             Krippendorpf's Alpha. Can be NOMINAL, ORDINAL, or INTERVAL.
%   makeplot: Plot the data, 1 or 0 (default).

% Calculate alpha with hist3 approach (absolute values)
tic

% Check inputs
if nargin < 2
    error('Error: must have argument SCALE')
elseif nargin < 3
    makeplot = 0;
end

fprintf('This dataset has %i observers and %i observations.\n',size(dat, 1) , size(dat, 2) )
scale = lower(scale);

% Get variables
allvals = unique(dat(~isnan(dat)));
tdim = size(dat, 2);

Y = repmat(1:size(dat,2), size(dat,1), 1);
XY = [dat(:), Y(:)];
dO = hist3(XY, {allvals, 1:tdim});
dE = sum(dO,2);
n__ = sum(sum(dO));              %length(dat(:));
nu_ = sum(dO,1);                 %size(dat,1);

% Optional plot
if makeplot
    figure;
    subplot(1,20,1:3); plot(dE, allvals); axis tight; hold on
    subplot(1,20,5:20); imagesc(1:size(dat,2),allvals,dO); colorbar
    set(gca,'YDir','normal')
end

% nominator
Zu = zeros(1, tdim);
for tt = 1:tdim
    if ~(nu_(tt) > 1)
        continue
    end
    % Find values
    Znucnuk = zeros(1, length(allvals)-1);  % Initiate with zeroes.
    vidx = find(dO(:,tt));                   % Find only index with observed values, all else will be zero anyway
    
    for ii = 1:length(vidx)-1
        idx = vidx(ii);
        c = allvals(idx);
        kvals = allvals(vidx(ii+1:end));
        switch scale
            case 'nominal'
                deltas = ~(kvals==c);
            case 'ordinal'
                deltas = nan(length(kvals),1);
                for gg = 1:length(kvals)
                    deltas(gg) = (sum(dE(ii:ii+gg)) - (dE(ii)+dE(ii+gg))/2)^2;
                end                    
            case 'interval'
                deltas = (kvals-c).^2;
        end

        nuc = dO(idx,tt);
        nuk = dO(vidx(ii+1:end), tt);
        Znucnuk(idx) = sum(nuc*nuk.*deltas);    
    end   
    Zu(tt) = sum(Znucnuk./(nu_(tt)-1));
end

% denominator
Zncnk = zeros(1, length(allvals));
for ii = 1:length(allvals)-1
    c = allvals(ii);
    kvals = allvals(ii+1:end);
    switch scale
        case 'nominal'
            deltas = ~(kvals==c);
        case 'ordinal'
            deltas = nan(length(kvals),1);
            for gg = 1:length(kvals)
                deltas(gg) = (sum(dE(ii:ii+gg)) - (dE(ii)+dE(ii+gg))/2)^2;
            end                    
        case 'interval'
            deltas = (kvals-c).^2;
    end        
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
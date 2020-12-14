function kripAlpha(dat, scale)
% Use as:
%   alpha = kripAlpha(X, method)
% Where 
%   dat:    N observers x M observations. For time series M = t
%   scale:  The method for calculating the error (i.e. delta) for
%           Krippendorpf's Alpha. Can be ...

% Check inputs


% Calculate alpha with hist3 approach (absolute values)
% tic
allvals = unique(dat);
tdim = 1:size(dat,2);

res = length(allvals); % data resolution

X = dat(:);
Y = repmat(1:size(dat,2),size(dat,1),1);
Y = Y(:);

XY = [X, Y];    %Data vale, time bin

% ll = linspace(min(dat(:)), max(dat(:)), res);
ll = allvals;
dO2 = hist3(XY, {ll, 1:size(dat,2)});

dE2 = sum(dO2,2);
totals2 = sum(dO2,1);

% figure;
% subplot(1,20,1:3); plot(dE2,ll); axis tight; hold on
% % plot(pgsum,xi2)
% subplot(1,20,5:20); imagesc(1:size(dat,2),ll,dO2); colorbar
% set(gca,'YDir','normal')

% P_vals = interp1(grid1, dE, allVals, 'spline');

%
n__ = sum(sum(dO2));              %length(dat(:));
nu_ = sum(dO2,1);                 %size(dat,1);

for tt = 1:length(tdim)
%     u = gridx2(tt);
    
    for ii = 1:length(allvals)-1
        c = allvals(ii);
        kvals = allvals(ii+1:end);
        deltas = (kvals-c).^2;
    
        % nominator
        nuc = dO2(ii,tt);
        nuk = dO2(ii+1:end, tt);
        Znucnuk(ii) = sum(nuc*nuk.*deltas);

        % denominator
        n_c = dE2(ii);
        n_k = dE2(ii+1:end);
        Zncnk(ii) = sum(n_c*n_k.*deltas);
        
    end
    Zu(tt) = sum(Znucnuk./(nu_(tt)-1));
end
toc

Do = sum(Zu);
De = sum(Zncnk) / (n__-1);

alphaH(aa) = 1 - (Do/De);

end


disp('done')
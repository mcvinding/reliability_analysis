function [alphaP] = alphaprime(dat, res)
% Method based on approximation for faster computation

%%
if ~exist('res', 'var')
    res = 0.01;     % 1% bin resolution
elseif res <= 0 || res >= 0.5
    error('Resultion must be larger than 0 and below 0.5')
end

%% Break
N = round(1/res);               % Number of bins
allvals = unique(dat(:));       % All unique values

%% Find histogram probability function over time
% Get joint density function over time
gridx1 = linspace(min(allvals), max(allvals), N);   % Data axis. Evenly spaces mased on max/min of data
tdim = size(dat, 2);                             % "Time" axis (copy from real data)

% X = interp1(gridx1, gridx1, dat(:), 'spline');
Y = repmat(1:size(dat, 2), size(dat, 1), 1);
XY = [dat(:), Y(:)];
dO = hist3(XY, {gridx1, 1:tdim}) ;
dE = sum(dO, 2);

% figure;
figure;
subplot(1,20,1:3); plot(dE, gridx1); axis tight; hold on
subplot(1,20,5:20); imagesc(1:size(dat,2),gridx1,dO); colorbar
set(gca,'YDir','normal')

n__ = sum(sum(dO));              %length(dat(:));
nu_ = sum(dO, 1);                 %size(dat,1);

tic

%% Interplation method
for tt = 1:tdim
    disp(tt)
%     u = gridx2(tt);
    tic
    intpt = interp1(gridx1, dO(:,tt), allvals, 'nearest');
    
    % Clean up
    X = gridx1; V = dO(:,tt); Xq = allvals;
    Vq = interp1(X,V,Xq, 'pchip'); 
    
    % divide by zero? No there must be a minimum of one valid value!
    c = sum(Vq)/nu_(tt);
    pVqall = Vq/(n__*c);
    Vqall = Vq/c;

    pVall = V/n__;

    subplot(1,2,1); plot(X,V,'o-',Xq,Vqall,':.')
    subplot(1,2,2); plot(X,pVall,'o-',Xq,pVqall,':.')
    
    % find test
    ZnucnukT = zeros(1, length(allvals)-1);  % Initiate with zeroes.
    vidx = find(Vqall);                      % Find only index with observed values, all else will be zero anyway
      
    for ii = 1:length(vidx)  
        idx = vidx(ii);
        c = allvals(idx);
        kvals = allvals(idx+1:end);
        deltas = (kvals-c).^2;

%         nuc = intpt(ii);
%         nuk = intpt(ii+1:end);     
        nuc = Vqall(idx);
        nuk = Vqall(idx+1:end);
        ZnucnukT2(idx) = sum(nuc*nuk.*deltas);  
    end
    Zu(tt) = sum(Znucnuk./(nu_(tt)-1));
    toc
end

% denominator
dEint = interp1(gridx1, dE, allvals, 'nearest');
X = gridx1; V = dE; Xq = allvals;
Vq = interp1(X,V,Xq, 'pchip'); 
plot(X,V,'o-',Xq,Vq,':.')

  
for ii = 1:length(allvals)-1
    c = allvals(ii);            % Real value
    kvals = allvals(ii+1:end);
    deltas = (kvals-c).^2;
        
    n_c = dEint(ii);
    n_k = dEint(ii+1:end);
    Zncnk(ii) = sum(n_c*n_k.*deltas);

end

toc

Do = sum(Zu);
De = sum(Zncnk) / (n__-1);

alphaP = 1 - (Do/De);



%% Smooth method
N = 100
% Get joint density function over time
gridx1 = linspace(min(dat(:)),max(dat(:)), N);   % Data axis. Evenly spaces mased on max/min of data
tdim = t;                                      % Time axis (copy from real data)

X = interp1(gridx1, gridx1, dat(:), 'spline');
Y = repmat(1:size(dat,2),size(dat,1),1);
Y = Y(:);
XY = [dat(:) Y];
dO = hist3(XY, {gridx1, 1:size(dat,2)}) ;

dE = sum(dO,2);

res_pct = 0.01; % smooth half-width resolution in pct
rx = 1/res_pct;
nx = N/rx;          % Kernel width

SIGMA = [nx 3];
dOblur = imgaussfilt(dO, SIGMA, 'Padding', 'replicate');
n__s = sum(sum(dOblur));              %length(dat(:));
nu_s = sum(dOblur, 1);                 %size(dat,1);

dEblur = sum(dOblur,2);

% figure;
figure;
subplot(1,20,1:3); plot(dEblur,gridx1); axis tight; hold on
subplot(1,20,5:20); imagesc(1:size(dat,2),gridx1,dOblur); colorbar
set(gca,'YDir','normal')

%% Calculate
for tt = 1:length(tdim)
%     u = gridx2(tt);
    
%     xt = Iblur(:,tt);
   intpt = interp1(gridx1, dOblur(:,tt), allvals, 'pchip'); 
   
    X = gridx1; V = dOblur(:,tt); Xq = allvals;
    Vq = interp1(X,V,Xq, 'nearest'); 
    plot(X,V,'o-',Xq,Vq,':.')

    for ii = 1:length(allvals)-1
        c = allvals(ii);            % Real value
        kvals = allvals(ii+1:end);
        deltas = (kvals-c).^2;
        
        nuc = intpt(ii);
        nuk = intpt(ii+1:end);     
        Znucnuk(ii) = sum(nuc*nuk.*deltas);
        
    end
    Zu(tt) = sum(Znucnuk./(nu_s(tt)-1));
end

% denominator
dEint = interp1(gridx1, dE, allvals, 'nearest');
dEblurint = interp1(gridx1, dEblur, allvals, 'nearest');

X = gridx1; V = dE; Xq = allvals;
Vq = interp1(X,V,Xq, 'spline'); 
plot(X,V,'o-',Xq,Vq,':.')

  
for ii = 1:length(allvals)-1
    c = allvals(ii);            % Real value
    kvals = allvals(ii+1:end);
    deltas = (kvals-c).^2;
        
    n_c = dEblurint(ii);
    n_k = dEblurint(ii+1:end);
    Zncnk(ii) = sum(n_c*n_k.*deltas);
end
        
        
toc

Do = sum(Zu);
De = sum(Zncnk) / (n__s-1);

alphaP = 1 - (Do/De)








%%

A = accumarray([xr, tdim], 1, [N length(tdim)])


kern_res = [dt*1, range(dat(:))/100];
    
gridx1 = linspace(min(dat(:)),max(dat(:)), N);   % Data axis. Evenly spaces mased on max/min of data
tdim = t;   

[x1,x2] = meshgrid(gridx1, tdim);
x1 = x1(:);
x2 = x2(:);
xi = [x2 x1];

tmp = dat';
yz = tmp(:);
xz = repmat(t, 1, size(dat,1));
Xz = [xz', yz];
[dO, xi1, kn] =  ksdensity(Xz,xi,'Bandwidth', kern_res, 'Function', 'pdf');

% Get probabilities
po = dO/sum(dO);
pO = reshape(po,length(tdim),length(gridx1))';
dO = reshape(dO,length(tdim),length(gridx1))';


% Get marginal
[dE, xi2] =  ksdensity(dat(:),gridx1,'Bandwidth', kern_res(2),'Function', 'pdf');

pE = dE/(sum(dE));
% pE = sum(pO,2);
dE2 = sum(dO,2);
pE2 = dE2/sum(dE2);
totals = sum(pO,1);

% Plot density
figure; 
subplot(1,20,1:3); plot(dE,xi2); axis tight; hold on
plot(dE2,xi2)
subplot(1,20,5:20); imagesc(tdim,gridx1,dO); colorbar
set(gca,'YDir','normal');
title(['Dens. aa=', num2str(aa)])

% Plot probability
% figure; 
% subplot(1,20,1:3); plot(pE,xi2); axis tight; hold on
% plot(pE2,xi2)
% subplot(1,20,5:20); imagesc(tdim,gridx1,pO); colorbar
% set(gca,'YDir','normal');
% title(['Prob. aa=', num2str(aa)])

n__ = sum(sum(dO));              %length(dat(:));
nu_ = sum(dO,1);                 %size(dat,1);

%% Calculate alpha
% Density approach
tic
allvals = unique(dat);

% nominator
for tt = 1:length(tdim)
    u = tdim(tt);
    
    for ii = 1:length(allvals)-1
        c = allvals(ii);
        kvals = allvals(ii+1:end);
        deltas = (kvals-c).^2;
    
        nuc = interp2(tdim, gridx1, dO, u, c, 'linear');
        nuk = interp2(tdim, gridx1, dO, u, kvals, 'linear');
        Znucnuk(ii) = sum(nuc*nuk.*deltas);
    end
    Zu(tt) = sum(Znucnuk./(nu_(tt)));
    Zu2(tt) = sum(Znucnuk);
end

% denominator
for ii = 1:length(allvals)-1
    c = allvals(ii);
    kvals = allvals(ii+1:end);
    deltas = (kvals-c).^2;
%     n_c = interp1(gridx1, dE2, c, 'linear');
%     n_k = interp1(gridx1, dE2, kvals, 'linear');

    n_c = interp1(gridx1, dE, c, 'linear');
    n_k = interp1(gridx1, dE, kvals, 'linear');

    Zncnk(ii) = sum(n_c*n_k.*deltas);
end

Do = sum(Zu);
Do2 = sum(Zu2);

De = sum(Zncnk); % / (n__-1);
De2 = sum(Zncnk) / (n__-1);


alphaD(aa+1) = 1 - (Do/De);
alphaD2(aa+1) = 1 - (Do2/De2);
alphaD3(aa+1) = 1 - (Do2/De);

T = toc;
fprintf('Calculating Alphaprime took %.2f sec\n', T)
end

%% Probability approach
tic
nu_ = sum(pO,1);
allvals = unique(dat);
for tt = 1:length(tdim)
    u = tdim(tt);
    
    for ii = 1:length(allvals)-1
        c = allvals(ii);
        kvals = allvals(ii+1:end);
        deltas = (kvals-c).^2;
    
        % nominator
        nuc = interp2(tdim, gridx1, pO, u, c, 'linear');
        nuk = interp2(tdim, gridx1, pO, u, kvals, 'linear');
        Znucnuk(ii) = sum(nuc*nuk.*deltas);

        % denominator
        n_c = interp1(gridx1, pE2, c, 'linear');
        n_k = interp1(gridx1, pE2, kvals, 'linear');
        Zncnk(ii) = sum(n_c*n_k.*deltas);
        
    end
    Zu(tt) = sum(Znucnuk);
end

% for tt = 1:length(tdim)
%     Zu(tt) = sum(Znucnuk./(nu_(tt)));
% end
    
    
Do = sum(Zu);
De = sum(Zncnk);

alphaP(aa) = 1 - (Do/De);

% timer
T = toc;
fprintf('Calculating Alphaprime took %.2f sec\n', T)

% end

%% hist3 approach (absolute values)
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


% %%
% figure;
% subplot(1,3,1); scatter(alpha1,alphaP)
% subplot(1,3,2); scatter(alpha1,alphaX2)
% subplot(1,3,3);scatter(alphaP,alphaX2)



%END       
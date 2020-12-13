% Test agreement  between two time-series

% Generate two time-series
fs      = 200; % Sampling frequency (samples per second) 
dt      = 1/fs; % seconds per sample 
StopTime = 1; % seconds 
t = (dt:dt:StopTime); % seconds 
Freq = 5; % Sine wave frequency (hertz) 

for aa = 1:10

x = sin(2*pi*Freq*t)*10+randn(1,length(t))*aa;
y = sin(2*pi*Freq*t)*10+randn(1,length(t))*aa;

figure(1);
subplot(2,5,aa); plot(t, x); hold on
subplot(2,5,aa); plot(t, y);
title(aa)

dat = [x; y];   % N observers vs M samples

%% Test with kriAplha function
tic
alpha1(aa) = kriAlpha(dat, 'interval')
toc
%% New method
% TO DO:
% * Option to specify resolution
N = 100;        % 100 = Default in linspace

%% Find joint probability function over time (
% Clean up this part. Naming is inconsitent and there are redundant parts.
kern_res = [dt*1, range(dat(:))/100];
      
% Get joint density function over time
gridx1 = linspace(min(dat(:)),max(dat(:)), N);   % Data axis. Evenly spaces mased on max/min of data
tdim = t;                                      % Time axis (copy from real data)
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
figure; 
subplot(1,20,1:3); plot(pE,xi2); axis tight; hold on
plot(pE2,xi2)
subplot(1,20,5:20); imagesc(tdim,gridx1,pO); colorbar
set(gca,'YDir','normal');
title(['Prob. aa=', num2str(aa)])

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

    n_c = interp1(gridx1, dE2, c, 'linear');
    n_k = interp1(gridx1, dE2, kvals, 'linear');

    Zncnk(ii) = sum(n_c*n_k.*deltas);
end

Do = sum(Zu);
Do2 = sum(Zu2);

De = sum(Zncnk) % / (n__-1);
De2 = sum(Zncnk) / (n__-1);


alphaD(aa) = 1 - (Do/De);
alphaD2(aa) = 1 - (Do2/De2);

dt = toc;
fprintf('Calculating Alphaprime took %.2f sec\n', dt)

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

for tt = 1:length(tdim)
    Zu(tt) = sum(Znucnuk./(nu_(tt)));
end
    
    
Do = sum(Zu);
De = sum(Zncnk);

alphaP(aa) = 1 - (Do/De);

% timer
dt = toc;
fprintf('Calculating Alphaprime took %.2f sec\n', dt)

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
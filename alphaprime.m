% Test agreement  between two time-series

% Generate two time-series
fs      = 200; % Sampling frequency (samples per second) 
dt      = 1/fs; % seconds per sample 
StopTime = 1; % seconds 
t = (dt:dt:StopTime); % seconds 
Freq = 5; % Sine wave frequency (hertz) 

x = sin(2*pi*Freq*t)*10 %+randn(1,length(t))*3;
y = sin(2*pi*Freq*t)*10 %+randn(1,length(t))*3;

figure; hold on
plot(t, x);
plot(t, y)

dat = [x; y];   % N observers vs M samples


%% Test with kriAplha function
% x = 1:20;
% y = rand(1,20);

% dat = [x; y];   % N observers vs M samples
tic
alpha1 = kriAlpha(dat, 'interval')
toc
%% New method
% TO DO:
% * Option to specify resolution
N = 100;        % Default in linspace

%% Find joint probability function over time
% Clean up this part. Naming is inconsitent and there are redundant parts.
kern_res = [dt*1, range(dat(:))/20];
      
% Get joint density function over time
gridx1 = linspace(min(dat(:)),max(dat(:)), N);   % Data axis. Evenly spaces mased on max/min of data
gridx2 = t;                                      % Time axis (copy from real data)
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x2 x1];

tmp = dat';
yz = tmp(:);
xz = repmat(t, 1, size(dat,1));
Xz = [xz', yz];
[dO, xi1, kn] =  ksdensity(Xz,xi,'Bandwidth', kern_res, 'PlotFcn','contour');

% Get probabilities
po = dO/sum(dO);
pO = reshape(po,length(gridx2),length(gridx1))';

% dO = reshape(dO,length(gridx2),length(gridx1))';


% Get marginal
[dE, xi2] =  ksdensity(dat(:),gridx1,'Bandwidth', kern_res(2));

pE = dE/(sum(dE));
pgsum = sum(pO,2);

totals = sum(dO,1);

% Plot joint 
figure;
subplot(1,20,1:3); plot(pE,xi2); axis tight; hold on
plot(pgsum,xi2)
subplot(1,20,5:20); imagesc(gridx2,gridx1,pO); colorbar
set(gca,'YDir','normal')

% Calculate alpha
for tt = 1:length(gridx2)
    u = gridx2(tt)
    for ii = 1:length(gridx1)-1
        c = gridx1(ii);
        kvals = gridx1(ii+1:end);
        deltas = (kvals-c).^2;
    
        % nominator
        nuc = pO(ii,tt);
        nuk = pO(ii+1:end, tt);
        sumnucnuk(ii) = sum(nuc*nuk.*deltas');

        % denominator
        n_c = pE(ii);
        n_k = pE(ii+1:end);
        sumncnk(ii) = sum(n_c*n_k.*deltas);
        
    end
    sumu(tt) = sum(sumnucnuk);
end

Do = sum(sumu.*(totals-1));
De = sum(sumncnk);

alphaX = 1 - Do/De

%% Coinciedse matrix approach
% Reshape the datamatric M (for M > 2)
allVals=unique(dat(:));
allVals=allVals(isfinite(allVals));

% for M = 1:size(dat,1)
%     
%     
% end

% Calculate dOck
kern_res = [range(allVals)/50, range(allVals)/50];
grid1 = linspace(min(allVals),max(allVals), N);  % Use this as approximation or all values as is now )mioght be inefficient)?

[x1,x2] = meshgrid(grid1, grid1);
x1 = x1(:);
x2 = x2(:);
xi = [x2 x1];
Xz = [dat(1,:)', dat(2,:)'];

% Confusion density
[ox, ~, kn] =  ksdensity(Xz,xi, 'PlotFcn','contour','Bandwidth', kern_res);
% ksdensity(Xz,xi, 'PlotFcn','contour','Bandwidth', kern_res)
po = ox/sum(ox);
pO = reshape(po,length(grid1),length(grid1))';

% imagesc(grid1,grid1,dX)

% Get marginal
[e, ~] =  ksdensity(allVals,grid1,'Bandwidth', min(kern_res(2)));
[e2, ~] =  ksdensity(allVals,grid2,'Bandwidth', min(kern_res(2)));

de = e/(sum(e));
pgsum1 = sum(pO,1);
pgsum2 = sum(pO,2);

% Plot joint 
figure;
subplot(5,5,[1 6 11, 16]); plot(pgsum1,grid1); axis tight; set(gca,'YDir','reverse')
subplot(5,5,[2:5 7:10, 12:15, 17:20]); imagesc(grid1,gridx1,pO);
subplot(5,5,22:25); plot(grid1,pgsum2);axis tight;

tic
P_vals = interp1(grid1, dE, allVals, 'spline');
Desum = zeros(length(allVals)-ii,1);
Dosum = zeros(length(allVals)-ii,1);
for ii = 1:length(allVals)-1
    c = allVals(ii);

%     De = zeros(length(allVals)-i-1,1);
%     Do = zeros(length(allVals)-i-1,1);

    kvals = allVals(ii+1:end);
    delta = (kvals-c).^2;
    Pck = interp2(grid1, grid1, pO, c, kvals, 'spline');
    dOC(ii) = sum(Pck.*delta);
    dEC(ii) = sum((P_vals(ii).*P_vals(ii+1:end)).*delta);
end

do = sum(dOC)
de = sum(dEC)
toc

alpha2 = 1-sum(do)/sum(de);




%%
    P_k = 
    
    for j = ii+1:length(allVals)
        fprintf('i=%i, j=%i\n', ii, j)
        c = allVals(ii);
        k = allVals(j);
%         P_k = interp1(grid1, dE, k, 'spline');
        Pck = interp2(grid1, grid1, pO, c, kvals, 'spline');
        
%     Puk = interp2(gridx2, gridx1, dM, tu, k, 'spline');
    
        De(j) = (P_c+P_k)*(c-k).^2;
        Do(j) = Pck*(c-k).^2;

%         tmp = 
%         do=do+coinMatr(c,k)*(allVals(c)-allVals(k)).^2;
%         de=de+nc(c)*nc(k)*(allVals(c)-allVals(k)).^2;
    end
    Desum(ii) = sum(De);
    Dosum(ii) = sum(Do);
end
toc

nc=sum(pO,1);
n=sum(nc);

alpha = 1-sum(Dosum)/sum(Desum);


%% Loop over repetitions

delta = (x-y).^2;   % Error fun for interval data (Krippendorff, 2018)


for u = 1:size(dat,2)
    c = dat(1,u);
    k = dat(2,u);
    tu = t(u);
    du = delta(u);
    
    % Find density of observed values
    % c = 2
    % k = 5
    % u = 0.100

    P_c = interp1(gridx1, pE, c, 'spline');
    P_k = interp1(gridx1, pE, k, 'spline');

    Puc = interp2(gridx2, gridx1, pO, tu, c, 'spline');
    Puk = interp2(gridx2, gridx1, pO, tu, k, 'spline');

    Do(u) = (Puc+Puk)*delta(u);
    De(u) 

end


%%



    
       
       
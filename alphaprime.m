% Test agreement  between two time-series

% Generate two time-series
fs      = 200; % Sampling frequency (samples per second) 
dt      = 1/fs; % seconds per sample 
StopTime = 1; % seconds 
t = (0:dt:StopTime); % seconds 
Freq = 5; % Sine wave frequency (hertz) 

x = sin(2*pi*Freq*t)*10+randn(1,length(t))*3;
y = sin(2*pi*Freq*t)*10+randn(1,length(t))*3;

figure; hold on
plot(t, x);
plot(t, y)

dat = [x; y];   % N observers vs M samples


%% Test with kriAplha function
% x = 1:20;
% y = rand(1,20);

dat = [x; y];   % N observers vs M samples
tic
alpha = kriAlpha(dat, 'interval')
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
gridx2 = t;                                   % Time axis (copy from real data)
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x2 x1];

tmp = dat';
yz = tmp(:);
xz = repmat(t, 1, size(dat,1));
Xz = [xz', yz];
[f, xi1, kn] =  ksdensity(Xz,xi,'Bandwidth', kern_res, 'PlotFcn','contour');

% Get probabilities
pf = f/sum(f);
dM = reshape(pf,length(gridx2),length(gridx1))';

% Get marginal
[g, xi2] =  ksdensity(dat(:),gridx1,'Bandwidth', kern_res(2));

pg = g/(sum(g));
pgsum = sum(dM,2);

% Plot joint 
figure;
subplot(1,20,1:3); plot(pg,xi2); axis tight; hold on
plot(pgsum,xi2)
subplot(1,20,5:20); imagesc(gridx2,gridx1,dM); colorbar
set(gca,'YDir','normal')

%% Sandbox
% Loop over repetitions

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

    P_c = interp1(gridx1, pg, c, 'spline');
    P_k = interp1(gridx1, pg, k, 'spline');

    Puc = interp2(gridx2, gridx1, dM, tu, c, 'spline');
    Puk = interp2(gridx2, gridx1, dM, tu, k, 'spline');

    Do(u) = (Puc+Puk)*delta(u);
    De(u) 

end



    
       
       
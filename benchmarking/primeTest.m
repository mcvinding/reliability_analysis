% Bechsmarking of the alphaprime function with 10% 1% and 0.1% binning 
% compared to absolute alpha
%
% Test agreement between two time-series: Sine-waves with random noise.

%% Settings
outdir = '/home/mikkel/reliability_analysis/benchmarking/output';

%% Generate two time-series : settings
nrun        = 100;
StopTime    = 10;                % seconds 
fs          = 1000;              % Sampling frequency (samples per second) 
dt          = 1/fs;              % seconds per sample 
t           = (dt:dt:StopTime);  % seconds 
Freq        = 5;                 % Sine wave frequency (hertz) 
err         = 0:0.1:2;           % Noise scaling

%% Init.
alpha = nan(nrun,length(err));
alphaP1 = nan(nrun,length(err));
alphaP2 = nan(nrun,length(err));
alphaP3 = nan(nrun,length(err));

%% Test 
for ii = 1:nrun
    fprintf('Run %i of %i... ', ii, nrun) 
    for jj = 1:length(err)       

        % Generate two time-series
        x = sin(2*pi*Freq*t)+randn(1,length(t))*err(jj);
        y = sin(2*pi*Freq*t)+randn(1,length(t))*err(jj);
        dat = [x; y];   % N observers vs M samples

        % Absolute method
        alpha(ii,jj) = kripAlpha(dat, 'interval'); disp('done')

        % Approximation method
        alphaP1(ii,jj) = alphaprime(dat, 0.1); disp('done')
        alphaP2(ii,jj) = alphaprime(dat, 0.01); disp('done')
        alphaP3(ii,jj) = alphaprime(dat, 0.001); disp('done')
        
    end
    disp('done')
end

% Save output
fprintf('Saving output... ')
save(fullfile(outdir,'primetest'), 'alpha','alphaP1','alphaP2','alphaP3')
disp('DONE')

%% Summaries
ciA0 = prctile(alpha, [5, 95]) - mean(alpha);
ciP1 = prctile(alphaP1, [5, 95]) - mean(alphaP1);
ciP2 = prctile(alphaP2, [5, 95]) - mean(alphaP2);
ciP3 = prctile(alphaP3, [5, 95]) - mean(alphaP3);

%% Plot test result
lw = 1.5;                 % LineWidth
ms = 12;                % MarkerSize
f1 = figure; hold on
errorbar(err+0.006, mean(alphaP3), ciP3(1,:), ciP3(2,:), 'b.-', 'MarkerSize',ms, 'lineWidth',lw);
errorbar(err+0.002, mean(alphaP2), ciP2(1,:), ciP2(2,:), 'b .--', 'MarkerSize',ms, 'lineWidth',lw);
errorbar(err-0.002, mean(alphaP1), ciP1(1,:), ciP1(2,:), 'b.-.', 'MarkerSize',ms, 'lineWidth',lw);
errorbar(err-0.006, mean(alpha), ciA0(1,:), ciA0(2,:), 'r.-', 'MarkerSize',ms, 'lineWidth',lw);
set(gcf, 'Position', [500, 500, 700, 800])
legend(['Alpha', char(39), ' (0.1%)'], ... 
       ['Alpha', char(39), ' (1%)'], ...
       ['Alpha', char(39), ' (10%)'], ...
       'Alpha')
xlim([0.00 2.05]);
ax = gca();
ax.LineWidth = lw;
xlabel('Noise factor'); ylabel('Alpha')
exportgraphics(f1, fullfile(outdir, 'primeComparison.jpg'), 'Resolution', 600); 
%close 

% %% Plot timing
% figure; hold on
% errorbar(1:10, mean(dT2), std(dT2)/2, 'bo-');
% errorbar(1:10, mean(dTp1), std(dTp1)/2, 'ro-');
% errorbar(1:10, mean(dTp2), std(dTp2)/2, 'ro--');
% errorbar(1:10, mean(dTp3), std(dTp3)/2, 'ro-.');
% legend('Alpha', ['Alpha', char(39), ' (10%)'], ['Alpha', char(39), ' (1%)'], ['Alpha', char(39), ' (0.1%)'])

%% Compare

kripAlpha([mean(alphaP3); mean(alpha)], 'interval')
kripAlpha([mean(alphaP2); mean(alpha)], 'interval')
kripAlpha([mean(alphaP1); mean(alpha)], 'interval')

% alphaprime([mean(alphaP1); mean(alpha)])
% alphaprime([mean(alphaP2); mean(alpha)])
% alphaprime([mean(alphaP3); mean(alpha)])

kripAlpha([alphaP3(:)'; alpha(:)'], 'interval')
kripAlpha([alphaP2(:)'; alpha(:)'], 'interval')
kripAlpha([alphaP1(:)'; alpha(:)'], 'interval')
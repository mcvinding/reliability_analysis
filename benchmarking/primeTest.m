% Benchmarking of the alphaprime function with 10% 1% and 0.1% binning 
% compared to absolute alpha
%
% Test agreement between two time-series: Sine-waves with random noise.

wrkdir = 'C:\Users\ncb623\reliability_analysis';
outdir = fullfile(wrkdir, '/benchmarking/output');
addpath(wrkdir)

%% Generate time-series settings
nrun        = 100;
StopTime    = 10;                % seconds 
fs          = 1000;              % Sampling frequency (samples per second) 
dt          = 1/fs;              % seconds per sample 
t           = (dt:dt:StopTime);  % seconds 
Freq        = 5;                 % Sine wave frequency (Hz) 
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
ciA0 = prctile(alpha, [2.5, 97.5]);
ciP1 = prctile(alphaP1, [2.5, 97.5]);
ciP2 = prctile(alphaP2, [2.5, 97.5]);
ciP3 = prctile(alphaP3, [2.5, 97.5]);

%% Plot test result
lw = 1;                 % LineWidth
ms = 4;                % MarkerSize

f1 = figure; hold on
errorbar(err-0.01, median(alpha), ciA0(1,:)-median(alpha), ciA0(2,:)-median(alpha), ...
    'ro--', 'MarkerSize', ms, 'lineWidth', lw);
errorbar(err-0.005, median(alphaP3), ciP3(1,:)-median(alphaP3), ciP3(2,:)-median(alphaP3), ...
    'o--', 'MarkerSize', ms, 'lineWidth', lw, 'color',"#0072BD");
errorbar(err+0.005, median(alphaP2), ciP2(1,:)-median(alphaP2), ciP2(2,:)-median(alphaP2), ...
    'o--', 'MarkerSize', ms, 'lineWidth', lw, 'color',"#D95319");
errorbar(err+0.01, median(alphaP1), ciP1(1,:)-median(alphaP1), ciP1(2,:)-median(alphaP1), ...
    'o--', 'MarkerSize', ms, 'lineWidth', lw, 'color',"#7E2F8E");

set(gcf, 'Position', [600, 400, 700, 600])
legend('Alpha', ...
       ['Alpha', char(39), ' (0.1%)'], ... 
       ['Alpha', char(39), ' (1%)'], ...
       ['Alpha', char(39), ' (10%)'])
xlim([-0.05 2.05]);
ax = gca();
ax.LineWidth = 1.5;
xlabel('Noise factor'); ylabel('Alpha')

exportgraphics(f1, fullfile(outdir, 'primeComparison.jpg'), 'Resolution', 600); 
%close 

%% Compare

kripAlpha([mean(alphaP3); mean(alpha)], 'ratio')
kripAlpha([mean(alphaP2); mean(alpha)], 'ratio')
kripAlpha([mean(alphaP1); mean(alpha)], 'ratio')

kripAlpha([alphaP3(:)'; alpha(:)'], 'interval')
kripAlpha([alphaP2(:)'; alpha(:)'], 'interval')
kripAlpha([alphaP1(:)'; alpha(:)'], 'interval')

for ii = 2:length(err)
    [~, pval1(ii)] = ttest(alphaP1(:,ii), alpha(:,ii), 'tail', 'both');
    [~, pval2(ii)] = ttest(alphaP2(:,ii), alpha(:,ii), 'tail', 'both');
    [~, pval3(ii)] = ttest(alphaP3(:,ii), alpha(:,ii), 'tail', 'both');

    aval1(ii) = kripAlpha([alphaP1(:,ii)'; alpha(:,ii)'], 'ratio');
    aval2(ii) = kripAlpha([alphaP2(:,ii)'; alpha(:,ii)'], 'ratio');
    aval3(ii) = kripAlpha([alphaP3(:,ii)'; alpha(:,ii)'], 'ratio');
end

%% %-deviation from real alpha
pctP1 = (alphaP1-alpha)./alpha*100;
mean(pctP1)
ciPctP1 = prctile(pctP1, [2.5, 97.5]);
pctP2 = (alphaP2-alpha)./alpha*100;
mean(pctP2)
ciPctP2 = prctile(pctP2, [2.5, 97.5]);
pctP3 = (alphaP3-alpha)./alpha*100;
mean(pctP3)
ciPctP3 = prctile(pctP3, [2.5, 97.5]);

%% Plot %
lw = 1.5;                 % LineWidth
ms = 4;         % MarkerSize
f2 = figure(); hold on
set(gcf, 'Position', [1300, 600, 700, 450])

subplot(2,1,2); hold on
errorbar(err-0.005, median(pctP3), ciPctP3(1,:)-median(pctP3), ciPctP3(2,:)-median(pctP3), ...
    'o--', 'MarkerSize',ms, 'lineWidth',1, 'color',"#0072BD");
errorbar(err+0.005, median(pctP2), ciPctP2(1,:)-median(pctP2), ciPctP2(2,:)-median(pctP2), ...
    'o--', 'MarkerSize',ms, 'lineWidth',1, 'color',"#D95319");
errorbar(err, median(pctP1), ciPctP1(1,:)-median(pctP1), ciPctP1(2,:)-median(pctP1), ...
    'o--', 'MarkerSize',ms, 'lineWidth',1, 'color',"#7E2F8E");
% set(gcf, 'Position', [1300, 600, 700, 200])
legend(['Alpha', char(39), ' (0.1%)'], ... 
       ['Alpha', char(39), ' (1%)'], ...
       ['Alpha', char(39), ' (10%)'], ...
       'Location', 'southwest')
xlim([0.05 2.05]); ylim([-12, 2])
ax = gca();
ax.LineWidth = lw;
xlabel('Noise factor'); 
ylabel('%-Error')

% exportgraphics(f2, fullfile(outdir, 'primeComparisonPct.jpg'), 'Resolution', 600); 

%% Difference from real alpha
difP1 = alphaP1-alpha;
mean(difP1)
ciDifP1 = prctile(difP1, [2.5, 97.5]);
difP2 = alphaP2-alpha;
mean(difP2)
ciDifP2 = prctile(difP2, [2.5, 97.5]);
difP3 = alphaP3-alpha;
mean(difP3)
ciDifP3 = prctile(difP3, [2.5, 97.5]);

%% Plot %
% lw = 1.5;                 % LineWidth
% ms = 12;                % MarkerSize
% f3 = figure; hold on
subplot(2,1,1); hold on

errorbar(err-0.005, median(difP3), ciDifP3(1,:)-median(difP3), ciDifP3(2,:)-median(difP3), ...
    'o--', 'MarkerSize',ms, 'lineWidth',1, 'color',"#0072BD");
errorbar(err+0.005, median(difP2), ciDifP2(1,:)-median(difP2), ciDifP2(2,:)-median(difP2), ...
    'o--', 'MarkerSize',ms, 'lineWidth',1, 'color',"#D95319");
errorbar(err, median(difP1), ciDifP1(1,:)-median(difP1), ciDifP1(2,:)-median(difP1), ...
    'o--', 'MarkerSize',ms, 'lineWidth',1, 'color',"#7E2F8E");
% set(gcf, 'Position', [1300, 600, 700, 200])
% legend(['Alpha', char(39), ' (0.1%)'], ... 
%        ['Alpha', char(39), ' (1%)'], ...
%        ['Alpha', char(39), ' (10%)'], ...
%        'Location', 'southeast')
xlim([-0.05 2.05]); ylim([-0.03, 0.005])
ax = gca();
ax.LineWidth = lw;
% xlabel('Noise factor'); 
ylabel('Difference')

% exportgraphics(f3, fullfile(outdir, 'primeComparisondif.jpg'), 'Resolution', 600); 
exportgraphics(f2, fullfile(outdir, 'primeComparisonDifPct.jpg'), 'Resolution', 600); 

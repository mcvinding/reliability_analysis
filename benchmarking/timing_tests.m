% Timing tests...
% Compare the time to  calculate Alpha with the various implemented
% methods. Plot average across N runs.
addpath('/home/mikkel/reliability_analysis')

%% Settings
outdir = '/home/mikkel/reliability_analysis/benchmarking/output';

%% Generate time-series settings
nrun        = 100;
StopTimes   = [100, 250, 500, round(logspace(3,5,17), -2)]; % 2.^(1:15); %, 15, 20, 25, 30, 40, 50];
fs          = 1;            % Sampling frequency (samples per second) 
dt          = 1/fs;         % seconds per sample 
Freq        = 5;

%% Init.
alphaOrg = nan(nrun, length(StopTimes));
alphaN2f = nan(nrun, length(StopTimes));
alphaPrm = nan(nrun, length(StopTimes));
alphaMat = nan(nrun, length(StopTimes));
dTOrg = nan(nrun, length(StopTimes));
dTN2f = nan(nrun, length(StopTimes));
dTPrm = nan(nrun, length(StopTimes));
dTMat = nan(nrun, length(StopTimes));

%% Run sims
for ii = 43:nrun
    fprintf('Run %i of %i...\n', ii, nrun)
    for aa = 1:length(StopTimes)
        t = (dt:dt:StopTimes(aa));
        fprintf('N = %i ... ', length(t))

        % Generate two time-series
        x = sin(2*pi*Freq*t)*10+randn(1,length(t))*5;
        y = sin(2*pi*Freq*t)*10+randn(1,length(t))*5;
        dat = [x; y];   % N observers vs M samples

        % New implementation
        if ~(length(t) > 75000)  % Above this and my PC will crash
            tic
            alphaOrg(ii,aa) = kripAlpha(dat, 'interval'); disp('done')
            dTOrg(ii,aa) = toc;
        end
        
        % New fast implementation
        tic
        alphaN2f(ii,aa) = kripAlphaN2fast(dat); disp('done')
        dTN2f(ii,aa) = toc;

        % Approximation method
        tic
        alphaPrm(ii,aa) = alphaprime(dat); disp('done')
        dTPrm(ii,aa) = toc;

        % Old MATLAB file
        if ~(length(t) > 500)  % Above this and my PC will crash
            tic
            alphaMat(ii,aa) = kriAlpha(dat, 'interval'); disp('done')
            dTMat(ii,aa) = toc;     
        end 

    end
end
disp('DONE')

% Save output
fprintf('Saving output... ')
save(fullfile(outdir,'timetest'), 'alphaOrg','alphaN2f','alphaPrm','alphaMat', ...
                                  'dTOrg','dTN2f','dTPrm','dTMat')
disp('DONE')
%% Summaries
ciOrg = prctile(dTOrg, [5, 95]) - mean(dTOrg, 'omitnan');
ciN2f = prctile(dTN2f, [5, 95]) - mean(dTN2f, 'omitnan');
ciPrm = prctile(dTPrm, [5, 95]) - mean(dTPrm, 'omitnan');
ciMat = prctile(dTMat, [5, 95]) - mean(dTMat, 'omitnan');

%% Plot
lw = 1.5;                 % LineWidth
ms = 5;                % MarkerSize

f1 = figure; hold on
set(gcf, 'Position', [500, 500, 700, 800])

errorbar(StopTimes*fs, mean(dTOrg, 'omitnan'), ciOrg(1,:), ciOrg(2,:), 'ko--', 'MarkerSize',ms, 'lineWidth',lw);
errorbar(StopTimes*fs, mean(dTN2f, 'omitnan'), ciN2f(1,:), ciN2f(2,:), 'bo--', 'MarkerSize',ms, 'lineWidth',lw);
errorbar(StopTimes*fs, mean(dTPrm, 'omitnan'), ciPrm(1,:), ciPrm(2,:), 'ro--', 'MarkerSize',ms, 'lineWidth',lw);
errorbar(StopTimes*fs, mean(dTMat, 'omitnan'), ciMat(1,:), ciMat(2,:), 'mo--', 'MarkerSize',ms, 'lineWidth',lw);
yline(300, '--'); yline(60, '--'); yline(10, '--')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend('Alpha', 'Alpha (N=2, fast)', ['Alpha', char(39)], 'Old MATLAB', 'location','northeast')
xlabel('Data points per observation'); ylabel('Time (s)')
ax = gca();
ax.LineWidth = lw;

exportgraphics(f1, fullfile(outdir, 'timeComparison.jpg'), 'Resolution', 600); 

%END
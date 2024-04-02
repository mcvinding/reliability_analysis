% Timing tests...
% Compare the time to  calculate Alpha with the various implemented
% methods. Plot average across N runs.

%% Init.
% Sample a time series
nrun = 100;
StopTimes = 2.^(1:15); %, 15, 20, 25, 30, 40, 50];
% fs          = 1000;            % Sampling frequency (samples per second) 
dt          = 1; %1/fs;            % seconds per sample 
Freq        = 5;

alphaOrg = nan(nrun, length(StopTimes));
alphaN2f = nan(nrun, length(StopTimes));
alphaPrm = nan(nrun, length(StopTimes));
alphaMat = nan(nrun, length(StopTimes));
dTOrg = nan(nrun, length(StopTimes));
dTN2f = nan(nrun, length(StopTimes));
dTPrm = nan(nrun, length(StopTimes));
dTMat = nan(nrun, length(StopTimes));

%% Run sims
for ii = 1:nrun
    fprintf('Run %i of %i...\n', ii, nrun)
    for aa = 1:length(StopTimes)
        t = (dt:dt:StopTimes(aa));
        fprintf('N = %i ... ', length(t))

        % Generate two time-series
        x = sin(2*pi*Freq*t)*10+randn(1,length(t))*5;
        y = sin(2*pi*Freq*t)*10+randn(1,length(t))*5;
        dat = [x; y];   % N observers vs M samples

%         % Test with kriAplha function
%         tic
%         alpha1b(ii,aa) = kriAlpha(dat, 'interval'); disp('done')
%         dT1b(ii,aa) = toc;

%         % New implementation
%         if ~(length(t) > 35000)  % Above this and my PC will crash
%             tic
%             alphaOrg(ii,aa) = kripAlpha(dat, 'interval'); disp('done')
%             dTOrg(ii,aa) = toc;
%         end
%         
%         % New fast implementation
%         tic
%         alphaN2f(ii,aa) = kripAlphaN2fast(dat); disp('done')
%         dTN2f(ii,aa) = toc;
% 
%         % Approximation method
%         tic
%         alphaPrm(ii,aa) = alphaprime(dat); disp('done')
%         dTPrm(ii,aa) = toc;

        % Old MATLAB file
        if ~(length(t) > 10000)  % Above this and my PC will crash
            tic
            alphaMat(ii,aa) = kriAlpha(dat, 'interval'); disp('done')
            dTMat(ii,aa) = toc;     
        end 

    end
end
disp('DONE')

%% Summaries

ciOrg = prctile(dTOrg, [5, 95]);
ciN2f = prctile(dTN2f, [5, 95]);
ciPrm = prctile(dTPrm, [5, 95]);

%% Plot

figure; hold on
errorbar(StopTimes*fs, mean(dTOrg), ciOrg(1,:), ciOrg(2,:), 'ko-');
errorbar(StopTimes*fs, mean(dTN2f), ciN2f(1,:), ciN2f(2,:), 'bo-');
errorbar(StopTimes*fs, mean(dTPrm), ciPrm(1,:), ciPrm(2,:), 'ro-');
% set(gca, 'YScale', 'log')
set(gcf, 'Position', [500, 500, 700, 500])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend('Alpha', 'Alpha (N=2, fast)', ['Alpha', char(39)], 'location','northwest')
xlabel('Data points per observation'); ylabel('Time (s)')
print('comparion_time','-dpng')
% Test agreement  between two time-series: Sine-waves with random noise

%% Generate two time-series : settings
fs          = 200;               % Sampling frequency (samples per second) 
dt          = 1/fs;              % seconds per sample 
StopTime    = 1;                 % seconds 
t           = (dt:dt:StopTime);  % seconds 
Freq        = 5;                 % Sine wave frequency (hertz) 

%% Test 
nrun = 100;

for ii = 1:nrun
    fprintf('Run %i of %i... ', ii, nrun) 
    for aa = 1:10
        % Generate two time-series
        x = sin(2*pi*Freq*t)*10+randn(1,length(t))*aa;
        y = sin(2*pi*Freq*t)*10+randn(1,length(t))*aa;
        dat = [x; y];   % N observers vs M samples

%         % Plot data
%         % figure(1);
%         subplot(2,5,aa); 
%         plot(t, x); hold on
%         % subplot(2,5,aa); 
%         plot(t, y);
%         title(aa)

%         % Test with kriAplha function (https://se.mathworks.com/matlabcentral/fileexchange/36016-krippendorff-s-alpha)
%         tic
%         alpha1(ii,aa) = kriAlpha(dat, 'interval'); disp('done')
%         dT1(ii,aa) = toc;

        % New implementation
        tic
        alpha2(ii,aa) = kripAlpha(dat, 'interval'); disp('done')
        dT2(ii,aa) = toc;

        % Approximation method
        tic
        alphap1(ii,aa) = alphaprime(dat, 0.01); disp('done')
        dTp1(ii,aa) = toc;
        
        tic
        alphap2(ii,aa) = alphaprime(dat, 0.1); disp('done')
        dTp2(ii,aa) = toc;
        
        tic
        alphap3(ii,aa) = alphaprime(dat, 0.001); disp('done')
        dTp3(ii,aa) = toc;
        
    end
    disp('done')
end
disp('DONE')

%% Plot
% Plot test result
figure; hold on
errorbar((0.1:0.1:1.0)-0.0075, mean(alpha2), std(alpha2)/2, 'r.-', 'MarkerSize',12);
errorbar((0.1:0.1:1.0)-0.0025, mean(alphap1), std(alphap1)/2, 'b.-', 'MarkerSize',12);
errorbar((0.1:0.1:1.0)+0.0025, mean(alphap2), std(alphap1)/2, 'b.--', 'MarkerSize',12);
errorbar((0.1:0.1:1.0)+0.0075, mean(alphap3), std(alphap1)/2, 'b.-.', 'MarkerSize',12);
set(gcf, 'Position', [500, 500, 700, 500])
legend('Alpha', ['Alpha', char(39), ' (1%)'], ['Alpha', char(39), ' (10%)'], ['Alpha', char(39), ' (0.1%)'])
xlim([0.05 1.05]);
xlabel('Noise factor'); ylabel('Alpha')
print('comparion1','-dpng')

% Plot timing
figure; hold on
errorbar(1:10, mean(dT2), std(dT2)/2, 'bo-', 'filled');
errorbar(1:10, mean(dTp1), std(dTp)/2, 'ro-');
errorbar(1:10, mean(dTp2), std(dTp)/2, 'ro--');
errorbar(1:10, mean(dTp3), std(dTp)/2, 'ro-.');


%% Time test
nrun = 100;
StopTimes = [1, 2, 4, 8, 10, 15, 20, 25, 30];
fs          = 1000;               % Sampling frequency (samples per second) 
dt          = 1/fs;              % seconds per sample 

for ii = 1:nrun
    fprintf('Run %i of %i...\n', aa, nrun)
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

        % New implementation
        tic
        alpha2b(ii,aa) = kripAlpha(dat, 'interval'); disp('done')
        dT2b(ii,aa) = toc;
        
        % New fast implementation
        tic
        alpha3b(ii,aa) = kripAlphaN2fast(dat); disp('done')
        dT3b(ii,aa) = toc;

        % Approximation method
        tic
        alphapb(ii,aa) = alphaprime(dat); disp('done')
        dTpb(ii,aa) = toc;
    end
end
disp('DONE')

%%
figure; hold on
errorbar(1:10, mean(dT2b), std(dT2b)/2, 'bo-');
errorbar(1:10, mean(dTpb), std(dTpb)/2, 'ro-')
errorbar(1:10, mean(dT3b), std(dTpb)/2, 'ko-')

set(gca, 'YScale', 'log')

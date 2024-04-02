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
        alphap1(ii,aa) = alphaprime(dat, 0.1); disp('done')
        dTp1(ii,aa) = toc;
        
        tic
        alphap2(ii,aa) = alphaprime(dat, 0.01); disp('done')
        dTp2(ii,aa) = toc;
        
        tic
        alphap3(ii,aa) = alphaprime(dat, 0.001); disp('done')
        dTp3(ii,aa) = toc;
        
    end
    disp('done')
end
disp('DONE')

%% Plot test result
figure; hold on
errorbar((0.1:0.1:1.0)-0.0030, mean(alpha2), std(alpha2)/2, 'r.-', 'MarkerSize',12);
errorbar((0.1:0.1:1.0)-0.0010, mean(alphap1), std(alphap1)/2, 'b.-', 'MarkerSize',12);
errorbar((0.1:0.1:1.0)+0.0010, mean(alphap2), std(alphap2)/2, 'b.--', 'MarkerSize',12);
errorbar((0.1:0.1:1.0)+0.0030, mean(alphap3), std(alphap3)/2, 'b.-.', 'MarkerSize',12);
set(gcf, 'Position', [500, 500, 700, 500])
legend('Alpha', ['Alpha', char(39), ' (10%)'], ['Alpha', char(39), ' (1%)'], ['Alpha', char(39), ' (0.1%)'])
xlim([0.05 1.05]);
xlabel('Noise factor'); ylabel('Alpha')
print('comparion1','-dpng'); 
%close 

%% Plot timing
figure; hold on
errorbar(1:10, mean(dT2), std(dT2)/2, 'bo-');
errorbar(1:10, mean(dTp1), std(dTp1)/2, 'ro-');
errorbar(1:10, mean(dTp2), std(dTp2)/2, 'ro--');
errorbar(1:10, mean(dTp3), std(dTp3)/2, 'ro-.');
legend('Alpha', ['Alpha', char(39), ' (10%)'], ['Alpha', char(39), ' (1%)'], ['Alpha', char(39), ' (0.1%)'])



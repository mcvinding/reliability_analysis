% Benchmarking errors
addpath('/home/mikkel/reliability_analysis')

%% Ordinal/binary data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1000;
Nperm = 100;
err = 0.0:0.005:0.5;  %  Percent disagreement

%% One-sample case
truedat = [ones(1,N*0.2), zeros(1,N*0.82)];
alph_ord1 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        errdat = truedat;
        rx = randperm(N);
        slct = rx(1:round(N*err(ee)));
    
        for ii = 1:length(slct)
            if errdat(slct(ii)) == 1
               errdat(slct(ii)) = 0;
            elseif errdat(slct(ii)) == 0
               errdat(slct(ii)) = 1;
            end
        end
    
        dat = [truedat; errdat];
        
        alph_ord1(ss,ee) = kripAlpha(dat, 'ordinal');
    end
end

% figure(1)
% scatter(err, alph1, 'b'); hold on
% plot(err, alph1, 'b')

figure(1); hold on
errorbar(err, mean(alph_ord1), std(alph_ord1), 'b', 'MarkerSize',12);
yline()

set(gcf, 'Position', [500, 500, 700, 500])

%% Independent samples
alph_ord2 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        dat1 = truedat;
        dat2 = truedat;
    
        rx1 = randperm(N);
        rx2 = randperm(N);
        slct1 = rx1(1:round(N*err(ee)));
        slct2 = rx2(1:round(N*err(ee)));
    
        for ii = 1:length(slct1)
            if dat1(slct1(ii)) == 1
               dat1(slct1(ii)) = 0;
            elseif dat1(slct1(ii)) == 0
               dat1(slct1(ii)) = 1;
            end
        end
    
        for ii = 1:length(slct2)
            if dat2(slct2(ii)) == 1
               dat2(slct2(ii)) = 0;
            elseif dat2(slct2(ii)) == 0
               dat2(slct2(ii)) = 1;
            end
        end
        dat = [dat1; dat2];
        
        alph_ord2(ss, ee) = kripAlpha(dat, 'ordinal');
    end
end

% scatter(err, alph2, 'r'); hold on
% plot(err, alph2, 'r')

figure(1); hold on
errorbar(err, mean(alph_ord2), std(alph_ord2), 'r', 'MarkerSize',12);


%% Interval data

%Generate two time-series : settings
fs       = 1000;             % Sampling frequency (samples per second) 
dt       = 1/fs;             % seconds per sample 
StopTime = 1;                % seconds 
t        = (dt:dt:StopTime); % seconds 
Freq     = 5;                % Sine wave frequency (hertz) 
err      = 0:0.1:2;          % Noise scaling
Nperm    = 100;

truedat = sin(2*pi*Freq*t);

%% One-sample: interval + added noise

alph_int1 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        errdat = truedat + randn(1,length(t))*err(ee);
    
        dat = [truedat; errdat];
        
        alph_int1(ss,ee) = kripAlpha(dat, 'interval');
    end
end

figure(2); hold on
errorbar(err, mean(alph_int1), std(alph_int1), 'b', 'MarkerSize',12);

%% Two-sample: interval + added noise
alph_int2 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        dat1 = truedat + randn(1,length(t))*err(ee);
        dat2 = truedat + randn(1,length(t))*err(ee);

        dat = [dat1; dat2];
        
        alph_int2(ss,ee) = kripAlpha(dat, 'interval');
    end
end
disp('done')

figure(2); hold on
errorbar(err, mean(alph_int2), std(alph_int2), 'r', 'MarkerSize',12);

%% One-sample: interval + replacement noise (no scaling)
err = 0.0:0.01:0.5;  %  Percent disagreement

alph_int3 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        errdat = truedat;
        rx = randperm(length(errdat));
        slct = rx(1:round(length(errdat)*err(ee)));
        errdat(slct) = (rand(1, length(slct))-0.5)*2;
        dat = [truedat; errdat];
        
        alph_int3(ss,ee) = kripAlpha(dat, 'interval');
    end
end
disp('done')

% figure; plot(truedat); hold on; plot(errdat)

figure(3); hold on
errorbar(err, mean(alph_int3), std(alph_int3), 'b', 'MarkerSize',12);

%% two-sample: interval + replacement noise (no scaling)
alph_int4 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        dat1 = truedat;
        dat2 = truedat;
        rx1 = randperm(length(errdat));
        rx2 = randperm(length(errdat));

        slct1 = rx1(1:round(length(errdat)*err(ee)));
        slct2 = rx2(1:round(length(errdat)*err(ee)));


        dat1(slct1) = (rand(1, length(slct1))-0.5)*2;
        dat2(slct1) = (rand(1, length(slct2))-0.5)*2;

        dat = [dat1; dat2];
        
        alph_int4(ss,ee) = kripAlpha(dat, 'interval');
    end
end
disp('done')

figure(3); hold on
errorbar(err, mean(alph_int4), std(alph_int4), 'r--', 'MarkerSize',12);


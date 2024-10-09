% Benchmarking: test how values of alpha are realted to %-agreement for
% ordinal, 
wrkdir = 'C:\Users\ncb623\reliability_analysis';
outdir = fullfile(wrkdir, '/benchmarking/output');
addpath(wrkdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nominal/binary data
N = 1000;
Nperm = 100;
errPct = 0.0:0.01:0.5;  %  Percent disagreement

truedat = [ones(1,N*0.05), zeros(1,N*0.95)];

%% Nominal: One-sample
alph_nom1 = zeros(Nperm, length(errPct));

for ee = 1:length(errPct)
    for ss = 1:Nperm
        errdat = truedat;
        rx = randperm(N);
        slct = rx(1:round(N*errPct(ee)));
    
        for ii = 1:length(slct)
            if errdat(slct(ii)) == 1
               errdat(slct(ii)) = 0;
            elseif errdat(slct(ii)) == 0
               errdat(slct(ii)) = 1;
            end
        end
    
        X = [truedat; errdat];
        
        alph_nom1(ss,ee) = kripAlpha(X, 'nominal');
    end
end
disp('done')

% figure(1)
% scatter(err, alph1, 'b'); hold on
% plot(err, alph1, 'b')

%% Nominal: Independent samples
alph_nom2 = zeros(Nperm, length(errPct));

for ee = 1:length(errPct)
    for ss = 1:Nperm
        dat1 = truedat;
        dat2 = truedat;
    
        rx1 = randperm(N);
        rx2 = randperm(N);
        slct1 = rx1(1:round(N*errPct(ee)));
        slct2 = rx2(1:round(N*errPct(ee)));
    
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

        X = [dat1; dat2];
        
        alph_nom2(ss, ee) = kripAlpha(X, 'nominal');
    end
end
disp('done')

% scatter(err, alph2, 'r'); hold on
% plot(err, alph2, 'r')

%% Plot
figure(1); hold on
errorbar(errPct, mean(alph_nom1), std(alph_nom1), 'b-', 'MarkerSize',12);
errorbar(errPct, mean(alph_nom2), std(alph_nom2), 'r-', 'MarkerSize',12);
set(gcf, 'Position', [500, 500, 700, 500])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ordinal data
N = 1000;
Nperm = 100;
errPct = 0.0:0.02:1;  %  Percent disagreement
ordVals = [1,2,3,4,5,6,7,8,9,10];

truedat = repmat(ordVals, 1, N/length(ordVals));
truedat = truedat(randperm(N));

%% Ordinal: One-sample
alph_ord1 = zeros(Nperm, length(errPct));
for ee = 1:length(errPct)
    for ss = 1:Nperm
        errdat = truedat;
        rx = randperm(N);
        slct = rx(1:round(N*errPct(ee)));

        for ii = 1:length(slct)
            tmp = ordVals(~(ordVals==errdat(slct(ii))));
            tmp = tmp(randperm(length(tmp)));
            errdat(slct(ii)) = tmp(1);
        end

        X = [truedat; errdat];
        
        alph_ord1(ss,ee) = kripAlpha(X, 'ordinal');
    end
end
disp('done')

%% Ordinal: independent sample
alph_ord2 = zeros(Nperm, length(errPct));
for ee = 1:length(errPct)
    for ss = 1:Nperm
        dat1 = truedat;
        dat2 = truedat;

        rx1 = randperm(N);
        rx2 = randperm(N);

        slct1 = rx1(1:round(N*errPct(ee)));
        slct2 = rx2(1:round(N*errPct(ee)));

        for ii = 1:length(slct1)
            tmp = ordVals(~(ordVals==dat1(slct1(ii))));
            tmp = tmp(randperm(length(tmp)));
            dat1(slct1(ii)) = tmp(1);
        end

        for ii = 1:length(slct2)
            tmp = ordVals(~(ordVals==dat2(slct2(ii))));
            tmp = tmp(randperm(length(tmp)));
            dat2(slct2(ii)) = tmp(1);
        end

        X = [dat1; dat2];
        
        alph_ord2(ss,ee) = kripAlpha(X, 'ordinal');
    end
end
disp('done')

%% Plot

figure(2); hold on
errorbar(errPct, mean(alph_ord1), std(alph_ord1), 'b-', 'MarkerSize',12);
errorbar(errPct, mean(alph_ord2), std(alph_ord2), 'r-', 'MarkerSize',12);
set(gcf, 'Position', [500, 500, 700, 500])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interval data
%Generate two time-series : settings
fs       = 1000;             % Sampling frequency (samples per second) 
dt       = 1/fs;             % seconds per sample 
StopTime = 1;                % seconds 
t        = (dt:dt:StopTime); % seconds 
Freq     = 5;                % Sine wave frequency (hertz) 
err   = 0:0.1:2;             % Noise scaling for noise added

Nperm    = 100;

truedat = sin(2*pi*Freq*t);

%% One-sample: interval + added noise

alph_int1 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        errdatAdd = truedat + randn(1,length(t))*err(ee);
        errdatRep = truedat + randn(1,length(t))*err(ee);

        X = [truedat; errdat];
        
        alph_int1(ss,ee) = kripAlpha(X, 'interval');
    end
end
disp('done')


%% Two-sample: interval + added noise
alph_int2 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        dat1 = truedat + randn(1,length(t))*err(ee);
        dat2 = truedat + randn(1,length(t))*err(ee);

        X = [dat1; dat2];
        
        alph_int2(ss,ee) = kripAlpha(X, 'interval');
    end
end
disp('done')

%% Plot
figure(3); hold on
errorbar(err, mean(alph_int1), std(alph_int1), 'b', 'MarkerSize',12);
errorbar(err, mean(alph_int2), std(alph_int2), 'r', 'MarkerSize',12);
set(gcf, 'Position', [500, 500, 700, 500])

%% One-sample: interval + replacement noise (no scaling)
err = 0.0:0.01:0.5;  %  Percent disagreement

alph_int3 = zeros(Nperm, length(err));

for ee = 1:length(err)
    for ss = 1:Nperm
        errdat = truedat;
        rx = randperm(length(errdat));
        slct = rx(1:round(length(errdat)*err(ee)));
        errdat(slct) = (rand(1, length(slct))-0.5)*2;
        X = [truedat; errdat];
        
        alph_int3(ss,ee) = kripAlpha(X, 'interval');
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

        X = [dat1; dat2];
        
        alph_int4(ss,ee) = kripAlpha(X, 'interval');
    end
end
disp('done')

figure(3); hold on
errorbar(err, mean(alph_int4), std(alph_int4), 'r--', 'MarkerSize',12);


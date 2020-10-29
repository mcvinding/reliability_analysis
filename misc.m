%% Test for similarity 

%% Data from tutorial
datOrg = [ 71, 73, 86, 59, 71;   % weight ratings
          74, 80,101, 62, 83;
          76, 80, 93, 66, 77]';
datOrg(:,:,2)  = [166,160,187,161,172;  ...       % height rating
               171,170,174,1
               63,182; ...
               171,165,185,162,181]';

% target x rater x varibales/(features)
% 
datOrg = [ 71;   % weight ratings
          74;]';
datOrg(:,:,2)  = [166;       % height rating
                  171]';
datOrg(:,:,3)  = [42;       % height rating
                  40]';              
%            

%% Replicate irr function from R
% Input
n = size(datOrg,1);              % Number of targets (ns)
c = size(datOrg,3);              % Number of variables (nv)
b = size(datOrg,2);              % Number of judges (nr)

r = size(nchoosek(1:b, 2), 1);   % Number of "judge" combinations

doSS = 0; deSS = 0;
for ii = 1:c
    tmp = datOrg(:,:,ii);
    SSw = sum(var(tmp,[],2)/n)*n*(b-1);             % SS within targets
    SSt = var(reshape(tmp,1,numel(tmp)))*(n*b-1);   % 
    SSc = var(mean(tmp))*n*(b-1);
    
    doSS = doSS+SSw;
    deSS = deSS+((b-1)*SSt+SSc);
end

iota = 1-(b*doSS)/deSS;  %output

%% New data (n x c x b)
datNew = zeros(3,2,5);
datNew(:,1,:) = [ 71, 73, 86, 59, 71;   % weight ratings
                  74, 80,101, 62, 83;
                  76, 80, 93, 66, 77]
datNew(:,2,:) = [ 166,160,187,161,172;  % height rating
                  171,170,174,163,182;
                  171,165,185,162,181]

% Input
n = size(datNew,3);              % Number of targets (ns)
c = size(datNew,2);              % Number of variables (nv)
b = size(datNew,1);              % Number of "judges" (nr)

r = size(nchoosek(1:b, 2), 1);   % Number of "judge" combinations

doSS = 0; deSS = 0;
for ii = 1:c
    tmp = squeeze(datNew(:,ii,:));       % b x n
    SSw = sum(var(tmp)/n)*n*(b-1);
    SSt = var(reshape(tmp,1,numel(tmp)))*(n*b-1);
    SSc = var(mean(tmp,2))*n*(b-1);
    
    doSS = doSS+SSw;
    deSS = deSS+((b-1)*SSt+SSc);
end

iota = 1-(b*doSS)/deSS;  %output
fprintf('Iota = %.3f\n', iota)

%% Test 
x = 1:20;
y = 2:21;

dat = [x; y]

% Input
n = size(dat,3);              % Number of targets (ns)
c = size(dat,2);              % Number of variables (nv)
b = size(dat,1);              % Number of "judges" (nr)

r = size(nchoosek(1:b, 2), 1);   % Number of "judge" combinations

doSS = 0; deSS = 0;
for ii = 1:c
    tmp = squeeze(dat(:,ii,:));       % b x n
    SSw = sum(var(tmp)/n)*n*(b-1);
    SSt = var(reshape(tmp,1,numel(tmp)))*(n*b-1);
    SSc = var(mean(tmp,2))*n*(b-1);
    
    doSS = doSS+SSw;
    deSS = deSS+((b-1)*SSt+SSc);
end

iota = 1-(b*doSS)/deSS;  %output
fprintf('Iota = %.3f\n', iota)


%% Old
% observed (d0)

% for r < s
% Input
n = size(datOrg,1);              % Number of targets (ns)
c = size(datOrg,3);              % Number of variables (nv)
b = size(datOrg,2);              % Number of "judges" (nr)
comb = size(nchoosek(1:b, 2), 1);   % Number of "judge" combinations

ssrs = zeros(n,b-1,b);
for r = 1:b-1
    for s = (r+1):b
        for ii = 1:n
            x = squeeze(datOrg(ii,r,:));
            y = squeeze(datOrg(ii,s,:));
            ssrs(ii,r,s) = sum((x-y).^2);
        end
    end
end

d0 = sum(sum(sum(ssrs))) / (n*comb);
    
ssrs = zeros(r,s);
for r = 1:b-1
    for s = (r+1):b
        ssij = zeros(n);
        for ii = 1:n
            for jj = 1:n
                x = squeeze(datOrg(ii,r,:));
                y = squeeze(datOrg(jj,s,:));
                
                ssij(ii,jj) = sum((x-y).^2);
            end
        end    
        ssrs(r,s) = sum(sum(ssij));
    end
end

de = sum(sum(ssrs)) / ((n^2)*comb)

iota = 1 - (d0/de);
fprintf('Iota = %.3f\n', iota)


%%
euc = sqrt(sum((x-y).^2));
d0 = sum(euc.^2) / (n*r);

d0 = sum(sum((x-y).^2)) / (n*r);

% expected (de)
tmp1 = zeros(c);
for ii = 1:c
    for jj = 1:c
        tmp1(jj,ii) = sum((x(:,ii)-y(:,jj))^2);
    end
end

de = sum(sum(tmp1)) / (n^2 * r)

iota = 1 - (d0/de);
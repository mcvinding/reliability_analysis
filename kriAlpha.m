function alpha=kriAlpha(data,scale)
% alpha=kriAlpha(data,scale)
%     calculates Krippendorff's Alpha as a measure of inter-rater agreement
%     data: rate matrix, each row is a rater or coder, each column is a case
%     scale: level of measurement, supported are 'nominal', 'ordinal', 'interval'
%     missing values have to be coded as NaN or inf

% For details about Krippendorff's Alpha see:
% http://en.wikipedia.org/wiki/Krippendorff%27s_Alpha
% Hayes, Andrew F. & Krippendorff, Klaus (2007). Answering the call for a
%   standard reliability measure for coding data. Communication Methods and
%   Measures, 1, 77-89

% Copyright notice: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/36016/versions/1/previews/kriAlpha.m/index.html

if nargin~=2
    help kriAlpha
    error('Wrong number of input arguments.')
end

allVals=unique(data(:));
allVals=allVals(isfinite(allVals));

% coincidence matrix
coinMatr=nan(length(allVals));
for r=1:length(allVals)
    for c=r:length(allVals)
        val=0;
        fprintf('r=%i, c=%i\n', r, c)
        for d=1:size(data,2)
            %find number of pairs
            thisEx=data(:,d);
            thisEx=thisEx(isfinite(thisEx));
            numEntr=length(thisEx);
            numP=0;
            for p1=1:numEntr
                for p2=1:numEntr
                    if p1==p2
                        continue
                    end
                    if (thisEx(p1)==allVals(r) && thisEx(p2)==allVals(c))
                        numP=numP+1;
                    end
                end
            end
            if numP
                val=val+numP/(numEntr-1);
            end
        end
        coinMatr(r,c)=val;
        coinMatr(c,r)=val;
    end
end

nc=sum(coinMatr,2);
n=sum(nc);

% expected agreement
expMatr=nan(length(allVals));
for i=1:length(allVals)
    for j=1:length(allVals)
        if i==j
            val=nc(i)*(nc(j)-1)/(n-1);
        else
            val=nc(i)*nc(j)/(n-1);
        end
        expMatr(i,j)=val;
    end
end

% difference matrix
diffMatr=zeros(length(allVals));
for i=1:length(allVals)
    for j=i+1:length(allVals)
        if i~=j
            if strcmp(scale, 'nominal')
                val=1;
            elseif strcmp(scale, 'ordinal')
                val=sum(nc(i:j))-nc(i)/2-nc(j)/2;
                val=val.^2;
            elseif strcmp(scale, 'interval')
                val=(allVals(j)-allVals(i)).^2;
            else
                error('unknown scale: %s', scale);
            end
        else
            val=0;
        end
        diffMatr(i,j)=val;
        diffMatr(j,i)=val;
    end
end

% observed - expected agreement
do=0; de=0;
for c=1:length(allVals)
    for k=c+1:length(allVals)
        if strcmp(scale, 'nominal')
            do=do+coinMatr(c,k);
            de=de+nc(c)*nc(k);
        elseif strcmp(scale, 'ordinal')
            do=do+coinMatr(c,k)*diffMatr(c,k);
            de=de+nc(c)*nc(k)*diffMatr(c,k);
        elseif strcmp(scale, 'interval')
            do=do+coinMatr(c,k)*(allVals(c)-allVals(k)).^2;
            de=de+nc(c)*nc(k)*(allVals(c)-allVals(k)).^2;
        else
            error('unknown scale: %s', scale);
        end
    end
end
de=1/(n-1)*de;
alpha=1-do/de;
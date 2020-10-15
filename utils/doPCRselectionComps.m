function [outcomes] = doPCRselectionComps(params)

%alphaLevel = params.alphaLevel;
data = params.data;
responseNormal  = params.responseNormal;

if isfield(params, 'nOrthcomp')
    nOrthcomp = params.nOrthcomp;
else
    nOrthcomp = 1;
end

if isfield(params, 'ncomp')
    ncomp = params.ncomp;
else
    ncomp = 3; % OPTncomp= nOrthcomp + nPrediccomp
end


% nOrthcomp = 1; ncomp =8;
% alphaLevel = 0.05;
[Z,W,Pv,T] = dosc(data,responseNormal,nOrthcomp,1E-3);
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Z,responseNormal,ncomp);
fprintf('Cumulative ratios of variation that the predictive components captured.\n');
fprintf('#PredComp:\t'); 
for kk=1:size(PCTVAR,2)
    fprintf('\t%d',kk);     
end
fprintf('\n');
fprintf('X:\t'); 
tpCum = cumsum(PCTVAR(1,:),2);
for kk=1:length(tpCum)
    fprintf('\t%.4f',tpCum(kk));     
end
fprintf('\n');
tpCum = cumsum(PCTVAR(2,:),2);
fprintf('Y:\t');
for kk=1:length(tpCum)
    fprintf('\t%.4f',tpCum(kk));     
end
fprintf('\n');


outINFO = struct([]);
allScores = [data*W data*XL(:,1:ncomp)]; %from OPLS
for kk=1:size(allScores,2)
    tpLoadingID = 1:kk;
    [b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1)  allScores(:,tpLoadingID)]);
    %stats: R^2, F static, p-value, an estimate of the error variance
    tpNext = length(outINFO)+1;
    outINFO(tpNext).aic = size(data,1)*log(stats(4)) + 2*(1+size(tpLoadingID,2));
    outINFO(tpNext).bic = size(data,1)*log(stats(4)) + log(size(data,1))*(1+size(tpLoadingID,2));
    outINFO(tpNext).adjR2 = 1- (size(data,1)-1)/(size(data,1)-kk-1)*(1-stats(1));
    if kk > 1
        BtpLoadingID = 1:(kk-1);
        [Bb,Bbint,Br,Brint,Bstats] = regress(responseNormal,[ones(size(data,1),1)  allScores(:,BtpLoadingID)]);
        [fstat fpval] = computeFstatFR(sum(r.^2), sum(Br.^2), kk, kk-1, size(data,1));    
        if kk==2
            fprintf('ncomp = nOrthcomp (%d) + #PredComp\n', nOrthcomp);
        end
        fprintf('[ncomp %d] fstat:%.4f fpval:%.3f\n', kk, fstat, fpval);
    end    
    outINFO(tpNext).stats = stats;
    outINFO(tpNext).type = 'OPLS';
    outINFO(tpNext).loadingID = tpLoadingID;
end

[coeffs, scores, latent,tsquare] = pca(data);
score_data = scores;
loading_data = coeffs;

allScores = scores(:,1:6); %from PCA
for kk=1:size(allScores,2)
    tpLoadingID = 1:kk;
    [b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1)  allScores(:,tpLoadingID)]);
    tpNext = length(outINFO)+1;
    outINFO(tpNext).aic = size(data,1)*log(stats(4)) + 2*(1+size(tpLoadingID,2));
    outINFO(tpNext).bic = size(data,1)*log(stats(4)) + log(size(data,1))*(1+size(tpLoadingID,2));
    outINFO(tpNext).adjR2 = 1- (size(data,1)-1)/(size(data,1)-kk-1)*(1-stats(1));
    outINFO(tpNext).stats = stats;
    outINFO(tpNext).type = 'PCA';
    outINFO(tpNext).loadingID = tpLoadingID;
end

showOutcomes = false;
allAIC = zeros(length(outINFO),1);
allBIC = zeros(length(outINFO),1);
for pp=1:length(outINFO)
    if showOutcomes
        fprintf('aic:%.4f bic:%.4f adjR2:%.3f pval:%.4f %s %s\n', outINFO(pp).aic, outINFO(pp).bic, outINFO(tpNext).adjR2, outINFO(pp).stats(3), outINFO(pp).type, num2str(outINFO(pp).loadingID));
    end
    allAIC(pp,1) = outINFO(pp).aic;
    allBIC(pp,1) = outINFO(pp).bic;
end
[allAICsort IDaic] = sort(allAIC);
[allAICsort2 IDaicRank] = sort(IDaic);

[allBICsort IDbic] = sort(allBIC);
[allBICsort2 IDbicRank] = sort(IDbic);

%[allAIC IDaicRank]
%[allBIC IDbicRank]
%[IDaicRank IDbicRank]

% for pp=1:length(outINFO)
%     %fprintf('aic:%.4f bic:%.4f pval:%.4f %s %s\n', outINFO(pp).aic, outINFO(pp).bic, outINFO(pp).stats(3), outINFO(pp).type, num2str(outINFO(pp).loadingID));
%     %if ( outINFO(pp).stats(3) < alphaLevel)
%     fprintf('%d pval:%.4f %s %s\n', IDaicRank(pp) + IDbicRank(pp), outINFO(pp).stats(3), outINFO(pp).type, num2str(outINFO(pp).loadingID));
%     %end
% end

fileID = fopen('tpAICBIC.txt','w');
for pp=1:length(outINFO)
    fprintf(fileID, 'aic:%.4f bic:%.4f  adjR2:%.3f pval:%.4f %s %s\n', outINFO(pp).aic, outINFO(pp).bic, outINFO(tpNext).adjR2, outINFO(pp).stats(3), outINFO(pp).type, num2str(outINFO(pp).loadingID));
end
fclose(fileID);

% settle the optimal setting
OPTtype = 'OPLS';
idPCAll = [1 2 3];
OPTncomp = 3;

%OPTtype = 'PCA';
%idPCAll = [1 2];

outcomes.OPTtype = OPTtype;
outcomes.idPCAll = idPCAll;
outcomes.OPTncomp = OPTncomp;
outcomes.ncomp = OPTncomp; %for compatibility

function [fstat fpval] = computeFstatFR(SSEf, SSEr, doef, doer, n)
fstat = (SSEr-SSEf)/(doef-doer)/(SSEf/(n-doef-1));
fpval = (1-fcdf(fstat, doef-doer, n-doef-1));
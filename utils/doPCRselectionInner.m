function [outcomes] = doPCRselectionInner(params)

data = params.data;
responseNormal = params.responseNormal;
alphaLevel = params.alphaLevel;
alphaLevelP = params.alphaLevelP;
choosetype = params.choosetype;

if isfield(params, 'idPCAll')
    idPCAll = params.idPCAll; %for PCA. It will not be used later.
else
    idPCAll = [1 2]; %for PCA. It will not be used later.
end


if isfield(params, 'OPTtype')
    OPTtype = params.OPTtype;
else
    OPTtype = 'OPLS';
end

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
OPTncomp = ncomp;

if isfield(params, 'usePercentage')
    usePercentage = params.usePercentage;
else
    usePercentage = 0.30;
end

if isfield(params, 'meanIterForOne')
    meanIterForOne = params.meanIterForOne;
else
    meanIterForOne = 300;
end

if isfield(params, 'minITERset')
    minITERset = params.minITERset;
else
    minITERset = 100;    
end

if isfield(params, 'maxITERset')
    maxITERset = params.maxITERset;
else
    maxITERset = 8000;    
end


if isfield(params, 'OPLSoption')
    OPLSoption = params.OPLSoption;
else
    OPLSoption = 2;
end

if isfield(params, 'OPLSIncludeOrtho')
    OPLSIncludeOrtho = params.OPLSIncludeOrtho;
else
    OPLSIncludeOrtho = 1;
    %OPLSIncludeOrtho = 0;
end


if isfield(params, 'verbose')
    verbose = params.verbose;
else
    verbose = true;
    %OPLSIncludeOrtho = 0;
end

%do variable-wise regression
%NumSigRegr = 0;
RegPvalues = zeros(size(data,2),1);
RegCoeffs = zeros(size(data,2),1);
for kk=1:size(data,2)
    [tb,tbint,tr,trint,tstats] = regress(responseNormal,[ones(size(data,1),1)  data(:,kk)]);
    RegPvalues(kk) = tstats(3);
    RegCoeffs(kk) = tb(2);
end
if verbose
    fprintf('current choosetype:%s\n',params.choosetype);
    fprintf(' # initial variables with %.3f significance (nocorrection) :%d\n', alphaLevel, sum( RegPvalues < alphaLevel ));
    fprintf(' # initial variables with %.3e significance (bonferoni):%d\n', alphaLevel/size(data,1), sum( RegPvalues < alphaLevel/size(data,1) ) );
end

%use bonferoni for this time. It should be set by the policy
if strcmp(choosetype, 'nocorrection')
    idxAll = find(RegPvalues < alphaLevel);
else
    idxAll = find(RegPvalues < alphaLevel/size(data,1));
end


%alphaLevel = 0.05;
%alphaLevelP = 0.05; %for the permutation test
upPCRvar = zeros(size(data,2),1);
lwPCRvar = zeros(size(data,2),1);
nOrthcomp = 1;


MAXITER = max(min(round(1/usePercentage*meanIterForOne), maxITERset), minITERset);
OrthWVector = [];
%obtaining PCRegr. vector
switch OPTtype
    case 'OPLS'
        switch OPLSoption
            case 1 %: OSC-corrected PLS RegrVector
                [Z,W,Pv,T] = dosc(data,responseNormal,nOrthcomp,1E-3);
                OrthWVector = W;
                [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Z,responseNormal,OPTncomp);
                pcr_vector = BETA(2:end) - W*Pv'*BETA(2:end);
            case 2 %all: OSC->PLS->Reg
                [Z,W,Pv,T] = dosc(data,responseNormal,nOrthcomp,1E-3);
                OrthWVector = W;
                [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Z,responseNormal,OPTncomp);
                if OPLSIncludeOrtho
                    %%Previous version
                    %[b,bint,r,rint,statsR] = regress(responseNormal,[ones(size(data,1),1)  data*W data*XL(:,1:OPTncomp)]);
                    %pcr_vector = [W XL(:,1:OPTncomp)]*b(2:end);
                    [b,bint,r,rint,statsR] = regress(responseNormal,[ones(size(data,1),1)  data*W data*stats.W]);
                    pcr_vector = [W stats.W]*b(2:end);
                else
                    %[b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1)  data*XL(:,1:OPTncomp)]);
                    %pcr_vector = [XL(:,1:OPTncomp)]*b(2:end);
                    
                    %To check out BETA in plsgress, pcr_vector should equal
                    %BETA(2:end) as below
                    %[b,bint,r,rint,statsR] = regress(responseNormal,[ones(size(data,1),1)  Z*stats.W]);
                    %pcr_vector = [stats.W]*b(2:end);
                    %norm( pcr_vector - BETA(2:end) ) / norm(pcr_vector)
                    
                    [b,bint,r,rint,statsR] = regress(responseNormal,[ones(size(data,1),1)  data*stats.W]);
                    pcr_vector = [stats.W]*b(2:end);
                end
            case 3 %OSC-> PLS
                [Z,W,Pv,T] = dosc(data,responseNormal,nOrthcomp,1E-3);
                OrthWVector = W;
                [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Z,responseNormal,OPTncomp);
                pcr_vector = BETA(2:end);
            case 4 %OSC-> PCR
                [Z,W,Pv,T] = dosc(data,responseNormal,nOrthcomp,1E-3);
                OrthWVector = W;
                [tcoeffs, tscores, tlatent,ttsquare] = pca(Z);
                [b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1)  tscores(:,OPTncomp)]);
                pcr_vector = tcoeffs(:,1:OPTncomp)*b(2:end);
            case 5 %PLS
                [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(data,responseNormal,OPTncomp);
                pcr_vector = BETA(2:end);
            case 6 %PCR
                [tcoeffs, tscores, tlatent,ttsquare] = pca(data);
                [b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1)  tscores(:,OPTncomp)]);
                pcr_vector = tcoeffs(:,1:OPTncomp)*b(2:end);
        end
        
        %stats: R^2, F static, p-value, an estimate of the error variance
        
    case 'PCA'
        [coeffs, scores, latent,tsquare] = pca(data);
        score_data = scores;
        loading_data = coeffs;
        [b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1)  score_data(:,idPCAll)]);
        %stats: R^2, F static, p-value, an estimate of the error variance
        pcr_vector = loading_data(:, idPCAll)*b(2:end);
end


%tpAllPCRcoeff = zeros(MAXITER,1);
%tpAllPCRcoeff = zeros(MAXITER, length(idxAll));
if length(idxAll) > 1
    for ee=1:length(idxAll)
        tpImportPCRcoeff(ee).coeff = [];
        tpOrthWVector(ee).coeff = [];
    end
else
    tpImportPCRcoeff = [];
    tpOrthWVector = [];
end

ticID = tic;

for uu=1:MAXITER
    %CurVarNum = floor(unifrnd(0,size(data,2)))+1;
    %fprintf('[%d-th] Now doing variable %d\n', uu, CurVarNum);
    data_used = data;
    %tpMeanVarable = mean(data(:,CurVarNum));
    %tpStdVarable = std(data(:,CurVarNum));
    %data_used(:,CurVarNum) = normrnd(tpMeanVarable,tpStdVarable,[size(data,1),1]);
    
    randID = randperm(length(idxAll));
    tpIDmax = round(length(idxAll)*usePercentage);
    idUsedInIdxAll = randID(1:tpIDmax);
    randID_used = idxAll(idUsedInIdxAll);
    %randID_used(1) == idxAll(randID(1)) == idxAll(idUsedInIdxAll(1))
    for kk=1:length(randID_used)
        CurVarNum = randID_used(kk);
        data_used(:,CurVarNum) = data(randperm(size(data,1)),CurVarNum);
    end
    if verbose
        if mod(uu,2000) == 0
           fprintf('...continuing the %d-th iteration with MAXITER %d, OPLSoption %d \n',uu, MAXITER, OPLSoption);
        end
    end
    switch OPTtype
        case 'OPLS'
            switch OPLSoption
                case 1 %: OSC-corrected PLS RegrVector
                    [Z_used,W_used,Pv_used,T_used] = dosc(data_used,responseNormal,nOrthcomp,1E-3);
                    OrthWVector_used = W_used;
                    [XL_used,YL,XS,YS,BETA,PCTVAR,MSE,stats_used] = plsregress(Z_used,responseNormal,OPTncomp);
                    pcr_vector_used = BETA(2:end) - W_used*Pv_used'*BETA(2:end);

                case 2 %all: OSC->PLS->Reg
                    [Z_used,W_used,Pv_used,T_used] = dosc(data_used,responseNormal,nOrthcomp,1E-3);
                    OrthWVector_used = W_used;
                    [XL_used,YL,XS,YS,BETA,PCTVAR,MSE,stats_used] = plsregress(Z_used,responseNormal,OPTncomp);
                    
                    if OPLSIncludeOrtho
                        [b_used,bint_used,r_used,rint_used,stats_usedR] = regress(responseNormal,[ones(size(data,1),1)  data*W_used data*stats_used.W]);
                        pcr_vector_used = [W_used stats_used.W]*b_used(2:end);
                    else
                        [b_used,bint_used,r_used,rint_used,stats_usedR] = regress(responseNormal,[ones(size(data,1),1)  data*stats_used.W]);
                        pcr_vector_used = [stats_used.W]*b_used(2:end);
                    end                    
                
                case 3 %OSC-> PLS
                    [Z_used,W_used,Pv_used,T_used] = dosc(data_used,responseNormal,nOrthcomp,1E-3);
                    OrthWVector_used = W_used;
                    [XL_used,YL,XS,YS,BETA,PCTVAR,MSE,stats_used] = plsregress(Z_used,responseNormal,OPTncomp);
                    pcr_vector_used = BETA(2:end);
                case 4 %OSC-> PCR
                    [Z_used,W_used,Pv_used,T_used] = dosc(data_used,responseNormal,nOrthcomp,1E-3);
                    OrthWVector_used = W_used;
                    [tcoeffs_used, tscores_used, tlatent_used,ttsquare_used] = pca(Z_used);
                    [b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1)  tscores_used(:,OPTncomp)]);
                    pcr_vector_used = tcoeffs_used(:,1:OPTncomp)*b(2:end);
                case 5 %PLS
                    [XL_used,YL,XS,YS,BETA,PCTVAR,MSE,stats_used] = plsregress(data_used,responseNormal,OPTncomp);
                    pcr_vector_used = BETA(2:end);
                case 6 %PCR
                    [tcoeffs_used, tscores_used, tlatent_used,ttsquare_used] = pca(data_used);
                    [b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1)  tscores_used(:,OPTncomp)]);
                    pcr_vector_used = tcoeffs_used(:,1:OPTncomp)*b(2:end);
            end
        case 'PCA'
            [loading_used, scores_used, latent,tsquare] = pca(data_used);
            [b,bint,r,rint,stats] = regress(responseNormal,[ones(size(data,1),1) scores_used(:,idPCAll)]);
            pcr_vector_used = loading_used(:, idPCAll)*b(2:end);
    end
    
    doGraph = false;
    if doGraph
        figure(31); clf; hold on;
        plot(pcr_vector_used, 'r');
        plot(pcr_vector, 'b');
        plot(CurVarNum, 0, 'go');
    end
    
    for kk=1:length(idUsedInIdxAll)
        %randID_used(1) == idxAll(randID(1)) == idxAll(idUsedInIdxAll(1))
        CurVarNum = idxAll(idUsedInIdxAll(kk));
        %idUsedInIdxAll(kk) is the index in the array of idxAll
        tpImportPCRcoeff(idUsedInIdxAll(kk)).coeff = ...
            [tpImportPCRcoeff(idUsedInIdxAll(kk)).coeff; pcr_vector_used(CurVarNum)];
        %tpOrthWVector(idUsedInIdxAll(kk)).coeff = ...
        %    [tpOrthWVector(idUsedInIdxAll(kk)).coeff; OrthWVector_used(CurVarNum)];
    end
end
% end

elapsedTime = toc(ticID);
if verbose
    fprintf('...elapsed time:%.3f\n',elapsedTime);
end

outcomes.tpImportPCRcoeff = tpImportPCRcoeff;
outcomes.idxAll = idxAll;
outcomes.pcr_vector = pcr_vector;
outcomes.RegPvalues = RegPvalues;
outcomes.RegCoeffs = RegCoeffs;
outcomes.OrthWVector = OrthWVector;
outcomes.tpOrthWVector = tpOrthWVector;
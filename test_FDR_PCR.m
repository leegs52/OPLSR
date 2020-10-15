clc; clear;
load flodat.mat
Y = [];
X = [];

Y = [Y; temper30];
X = [X; spec30'];
Y = [Y; temper40];
X = [X; spec40'];
Y = [Y; temper50];
X = [X; spec50'];
Y = [Y; temper60];
X = [X; spec60'];
Y = [Y; temper70];
X = [X; spec70'];

% %downsampling
% randidx=randsample(110,80)
% X=X(randidx,:)
% Y=Y(randidx)

%MAXITER = 20000;
MAXITER = 200;
fdrTRsse = zeros(MAXITER,1);
fdrTRsst = zeros(MAXITER,1);
fdrTRR2 = zeros(MAXITER,1);
fdrLen = zeros(MAXITER,1);

fdrTEsse = zeros(MAXITER,1);
fdrTEsst = zeros(MAXITER,1);
fdrTER2 = zeros(MAXITER,1);

pcrTRsse = zeros(MAXITER,1);
pcrTRsst = zeros(MAXITER,1);
pcrTRR2 = zeros(MAXITER,1);
pcrLen = zeros(MAXITER,1);

pcrTEsse = zeros(MAXITER,1);
pcrTEsst = zeros(MAXITER,1);
pcrTER2 = zeros(MAXITER,1);

pcrTRsseO3 = zeros(MAXITER,1);
pcrTRsstO3 = zeros(MAXITER,1);
pcrTRR2O3 = zeros(MAXITER,1);
pcrLenO3 = zeros(MAXITER,1);

pcrTEsseO3 = zeros(MAXITER,1);
pcrTEsstO3 = zeros(MAXITER,1);
pcrTER2O3 = zeros(MAXITER,1);


%%
for uu=1:MAXITER
    if mod(uu,200) == 0
        fprintf('Now uu %d\n', uu);
    end
    trperc = 0.7;
    tpRandOrder = randperm(size(X,1));
    tpCutoff = round(size(X,1)*trperc);
    trID = tpRandOrder(1:tpCutoff);
    testID = tpRandOrder(tpCutoff+1:end);
    trX = X(trID,:);
    trY = Y(trID,:);
    
    teX = X(testID,:);
    teY = Y(testID,:);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alphaLevel = 0.01;
    %alphaLevel = 0.01;
    alphaLevelP = 0.05; %permutation    
    qlev = alphaLevel;
    ncomp = 3;
    nOrthcomp = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fdr
    data = trX;
    response = trY;
    %responseNormal=  trY;
    %nOrthcomp = 1;
    %ncomp = 8;
    %alphaLevel = 0.05;
    %alphaLevelP = 0.05;    
    %qlev = alphaLevel;
    %qlev = 0.2;
    
    chem_shift_FDR = (1:size(data,2))'; %p-by-1
    % use data' in p-by-n
    % response' in 1-by-n
    [sig_feature] = doFDRreg(data', response', qlev, chem_shift_FDR, 1);
    %[sig_feature] = doFDR(data', response', qlev, chem_shift_FDR, 1);
    
    %features found
    % test the model
    fdrselectid = sig_feature(:,end);

    tpfdrTRsse= 0;
    tpfdrTRsst = 0;
    tpfdrTRR2 = 0;
    tpfdrLen = 0;
    tpfdrTEsse = 0;
    tpfdrTEsst = 0;
    tpfdrTER2 = 0;
    
    if length(fdrselectid) > 0
        trXtruc = trX(:,fdrselectid);
        if size(trXtruc,1) > size(trXtruc,2)
            beta = inv(trXtruc'*trXtruc)*trXtruc'*trY;
        else
            %beta = pinv(trXtruc)*trY;
            beta = pinv(trXtruc'*trXtruc)*trXtruc'*trY;
        end
        trRes = trY - trXtruc*beta;
        tpfdrTRsse = sum(trRes.^2);
        tpfdrTRsst = sum(trY.^2);
        tpfdrTRR2 = 1 - tpfdrTRsse/tpfdrTRsst;
        tpfdrLen = length(fdrselectid);
        
        fprintf('[fdr] train R2:%.3f\n', tpfdrTRR2);
        
        teRes = teY - teX(:,fdrselectid)*beta;        
        tpfdrTEsse = sum(teRes.^2);
        tpfdrTEsst = sum(teY.^2);
        tpfdrTER2 = 1 - tpfdrTEsse/tpfdrTEsst;
        fprintf('[fdr] test R2:%.3f\n', tpfdrTER2);
    else
        fprintf('[fdr] No var found.\n');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pcr -- OPLSoption = 1
    data = trX;
    %nOrthcomp = 1;
    %ncomp = 3;
    %alphaLevel = 0.05;
    %alphaLevelP = 0.05;
    
    
    % step 1: selection of loadings starting with OPLS
    params.nOrthcomp = nOrthcomp;
    params.ncomp = ncomp;
    params.data = data;
    params.responseNormal = trY;
    
    OPTtype = 'OPLS';
    %idPCAll = [1 2 3];
    %OPTncomp = 3;
    OPLSoption = 1;
    
    % step 2: global permutation test; proceeding with permutation tests
    params.data = data;
    params.responseNormal = trY;
    params.alphaLevel = alphaLevel;
    params.alphaLevelP = alphaLevelP;
    %params.choosetype = 'bonferoni';
    params.choosetype = 'nocorrection';
    params.nOrthcomp = 1;
    %params.ncomp = OPTncomp;
    params.OPTtype = OPTtype;
    params.OPLSoption = OPLSoption;
    %params.idPCAll = idPCAll;
    params.minITERset = 30;
    params.maxITERset = 300;
    
    [outcomes2] = doPCRselectionInner(params);
    % step 3: summarize the result
    params.alphaLevelP = alphaLevelP;
    params.tpImportPCRcoeff = outcomes2.tpImportPCRcoeff;
    params.idxAll = outcomes2.idxAll;
    params.pcr_vector = outcomes2.pcr_vector;
    params.RegPvalues = outcomes2.RegPvalues;
    params.RegCoeffs = outcomes2.RegCoeffs;
    params.doSave = false;
    
    outcomes3 = doPCRselectionSummarize(params);
    % test the model
    pcr_vector = outcomes2.pcr_vector;
    pcrselectid = outcomes3.select;
    tppcrTRsse= 0;
    tppcrTRsst = 0;
    tppcrTRR2 = 0;
    tppcrLen = 0;
    tppcrTEsse = 0;
    tppcrTEsst = 0;
    tppcrTER2 = 0;
    if length(pcrselectid) >0
        trXtruc = trX(:,pcrselectid);
        if size(trXtruc,1) > size(trXtruc,2)
            beta = inv(trXtruc'*trXtruc)*trXtruc'*trY;
        else
            %beta = pinv(trXtruc)*trY;
            beta = pinv(trXtruc'*trXtruc)*trXtruc'*trY;
        end
        trRes = trY - trXtruc*beta;        
        %trRes = trY - trX(:,pcrselectid)*pcr_vector(pcrselectid);        
        
        tppcrTRsse = sum(trRes.^2);
        tppcrTRsst = sum(trY.^2);
        tppcrTRR2 = 1 - sum(trRes.^2)/sum(trY.^2);
        tppcrLen = length(pcrselectid);
        fprintf('[pcr] train R2:%.3f\n', tppcrTRR2);
        %fprintf('train R^2:%.3f\n', 1-sum(trRes.^2)/sum((trY).^2));
        
        %TeresponseNormal = zscore(teY);
        %teRes = teY - teX(:,pcrselectid)*pcr_vector(pcrselectid);
        teRes = teY - teX(:,pcrselectid)*beta;
        tppcrTEsse = sum(teRes.^2);
        tppcrTEsst = sum(teY.^2);
        tppcrTER2 = 1-sum(teRes.^2)/sum(teY.^2);
        fprintf('[pcr-%d] test R2:%.3f\n', OPLSoption, tppcrTER2);
        %fprintf('test R^2:%.3f\n', 1-sum(teRes.^2)/sum((teY).^2));
    else
        fprintf('[pcr-%d] No var found.\n',OPLSoption);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pcr --- OPLSoption = 3
    data = trX;
    %nOrthcomp = 1;
    %ncomp = 3;
    %alphaLevel = 0.05;
    %alphaLevelP = 0.05;
    
    
    % step 1: selection of loadings starting with OPLS
    params.nOrthcomp = nOrthcomp;
    params.ncomp = ncomp;
    params.data = data;
    params.responseNormal = trY;
    
    OPTtype = 'OPLS';
    %idPCAll = [1 2 3];
    %OPTncomp = 3;
    OPLSoption = 3;
    
    % step 2: global permutation test; proceeding with permutation tests
    params.data = data;
    params.responseNormal = trY;
    params.alphaLevel = alphaLevel;
    params.alphaLevelP = alphaLevelP;
    %params.choosetype = 'bonferoni';
    params.choosetype = 'nocorrection';
    params.nOrthcomp = 1;
    %params.ncomp = OPTncomp;
    params.OPTtype = OPTtype;
    params.OPLSoption = OPLSoption;
    %params.idPCAll = idPCAll;
    params.minITERset = 30;
    params.maxITERset = 300;
    
    [outcomes2] = doPCRselectionInner(params);
    % step 3: summarize the result
    params.alphaLevelP = alphaLevelP;
    params.tpImportPCRcoeff = outcomes2.tpImportPCRcoeff;
    params.idxAll = outcomes2.idxAll;
    params.pcr_vector = outcomes2.pcr_vector;
    params.RegPvalues = outcomes2.RegPvalues;
    params.RegCoeffs = outcomes2.RegCoeffs;
    params.doSave = false;
    
    outcomes3 = doPCRselectionSummarize(params);
    % test the model
    pcr_vector = outcomes2.pcr_vector;
    pcrselectidO3 = outcomes3.select;
    tppcrTRsseO3 = 0;
    tppcrTRsstO3 = 0;
    tppcrTRR2O3 = 0;
    tppcrLenO3 = 0;
    tppcrTEsseO3 = 0;
    tppcrTEsstO3 = 0;
    tppcrTER2O3 = 0;
    if length(pcrselectidO3) >0
        trXtruc = trX(:,pcrselectidO3);
        if size(trXtruc,1) > size(trXtruc,2)
            beta = inv(trXtruc'*trXtruc)*trXtruc'*trY;
        else
            %beta = pinv(trXtruc)*trY;
            beta = pinv(trXtruc'*trXtruc)*trXtruc'*trY;
        end
        trResO3 = trY - trXtruc*beta;        
        %trRes = trY - trX(:,pcrselectid)*pcr_vector(pcrselectid);        
        
        tppcrTRsseO3 = sum(trResO3.^2);
        tppcrTRsstO3 = sum(trY.^2);
        tppcrTRR2O3 = 1 - sum(trResO3.^2)/sum(trY.^2);
        tppcrLenO3 = length(pcrselectidO3);
        fprintf('[pcr-%d] train R2:%.3f\n', OPLSoption, tppcrTRR2O3);
        %fprintf('train R^2:%.3f\n', 1-sum(trRes.^2)/sum((trY).^2));
        
        %TeresponseNormal = zscore(teY);
        %teRes = teY - teX(:,pcrselectid)*pcr_vector(pcrselectid);
        teResO3 = teY - teX(:,pcrselectidO3)*beta;
        tppcrTEsseO3 = sum(teResO3.^2);
        tppcrTEsstO3 = sum(teY.^2);
        tppcrTER2O3 = 1-sum(teResO3.^2)/sum(teY.^2);
        fprintf('[pcr-%d] test R2:%.3f\n', OPLSoption, tppcrTER2O3);
        %fprintf('test R^2:%.3f\n', 1-sum(teRes.^2)/sum((teY).^2));
    else
        fprintf('[pcr-%d] No var found.\n', OPLSoption);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('FDR %d: tr %.3e \t te %.3f\n', length(fdrselectid), tpfdrTRR2, tpfdrTER2);
    fprintf('PCR-1 %d: tr %.3e \t te %.3f\n', length(pcrselectid), tppcrTRR2, tppcrTER2);
    fprintf('PCR-3 %d: tr %.3e \t te %.3f\n', length(pcrselectidO3), tppcrTRR2O3, tppcrTER2O3);
    
    fdrTRsse(uu,1) = tpfdrTRsse;
    fdrTRsst(uu,1) = tpfdrTRsst;
    fdrTRR2(uu,1) = tpfdrTRR2;
    fdrLen(uu,1) = tpfdrLen;
    
    fdrTEsse(uu,1) = tpfdrTEsse;
    fdrTEsst(uu,1) = tpfdrTEsst;
    fdrTER2(uu,1) = tpfdrTER2;
    
    pcrTRsse(uu,1) = tppcrTRsse;
    pcrTRsst(uu,1) = tppcrTRsst;
    pcrTRR2(uu,1) = tppcrTRR2;
    pcrLen(uu,1) = tppcrLen;
    
    pcrTEsse(uu,1) = tppcrTEsse;
    pcrTEsst(uu,1) = tppcrTEsst;
    pcrTER2(uu,1) = tppcrTER2;
    
    pcrTRsseO3(uu,1) = tppcrTRsseO3;
    pcrTRsstO3(uu,1) = tppcrTRsstO3;
    pcrTRR2O3(uu,1) = tppcrTRR2O3;
    pcrLenO3(uu,1) = tppcrLenO3;
    
    pcrTEsseO3(uu,1) = tppcrTEsseO3;
    pcrTEsstO3(uu,1) = tppcrTEsstO3;
    pcrTER2O3(uu,1) = tppcrTER2O3;
end

% save('NIRresult01.mat', ...
%     'fdrTRsse', 'fdrTRsst', 'fdrTRR2', 'fdrLen', 'fdrTEsse', 'fdrTEsst', 'fdrTER2', ...
%     'pcrTRsse', 'pcrTRsst', 'pcrTRR2', 'pcrLen', 'pcrTEsse', 'pcrTEsst', 'pcrTER2', ...
%     'pcrTRsseO3', 'pcrTRsstO3', 'pcrTRR2O3', 'pcrLenO3', 'pcrTEsseO3', 'pcrTEsstO3', 'pcrTER2O3', ...
%     'alphaLevel', 'qlev', 'ncomp', 'nOrthcomp');

%% 
load NIRresult01.mat
load NIRresult.mat
load NIRresult10.mat

%Q2
mean( fdrTER2(fdrLen~=0) )
mean( pcrTER2(pcrLen~=0) )
mean( pcrTER2O3(pcrLenO3~=0.5) )

%mse
mean( fdrTEsse(fdrLen~=0) )
mean( pcrTEsse(pcrLen~=0) )
mean( pcrTEsseO3(pcrLenO3~=0) )

%how many variables found
mean(fdrLen(fdrLen~=0))
mean(pcrLen(pcrLen~=0))
mean(pcrLenO3(pcrLenO3~=0))

%how many times variables are found
sum(fdrLen~=0)/length(fdrLen)
sum(pcrLen~=0)/length(pcrLen)
sum(pcrLenO3~=0)/length(pcrLenO3)

%% to save figures
[outcomes3.select  wl(outcomes3.select)]
%          381         960
%          382         961
%          383         962
%          384         963
%          385         964
%          386         965
%          387         966
%          388         967
%          389         968
%          390         969
%          470        1049
%          474        1053
%          476        1055
%          477        1056
%          481        1060
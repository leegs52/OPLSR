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

% randidx=randsample(110,80)
% X=X(randidx,:)
% Y=Y(randidx)

MAXITER = 20;

lassoTRsse = zeros(MAXITER,1);
lassoTRsst = zeros(MAXITER,1);
lassoTRR2 = zeros(MAXITER,1);
lassoLen = zeros(MAXITER,1);

lassoTEsse = zeros(MAXITER,1);
lassoTEsst = zeros(MAXITER,1);
lassoTER2 = zeros(MAXITER,1);

%%

% for uu=1:MAXITER
%     if mod(uu,20) == 0
%         fprintf('Now uu %d\n', uu);
%     end
%     trperc = 0.7;
%     tpRandOrder = randperm(size(X,1));
%     tpCutoff = round(size(X,1)*trperc);
%     trID = tpRandOrder(1:tpCutoff);
%     testID = tpRandOrder(tpCutoff+1:end);
%     trX = X(trID,:);
%     trY = Y(trID,:);
%     
%     teX = X(testID,:);
%     teY = Y(testID,:);
%     
%     data = trX;
%     response = trY;
%     k=0.1;
%     [sig_feature_r]=ridge(response, data, k, );
%     ridge_selectid = find(abs(sig_feature_r)>0.1);
%     
%     tpridgeTRsse= 0;
%     tpridgeTRsst = 0;
%     tpridgeTRR2 = 0;
%     tpridgeLen = 0;
%     tpridgeTEsse = 0;
%     tpridgeTEsst = 0;
%     tpridgeTER2 = 0;
%     
%     if length(ridge_selectid) > 0
%         trXtruc = trX(:,ridge_selectid);
%         beta=sig_feature_r(ridge_selectid);
%         
%         trRes = trY - (beta(1)+trXtruc*beta);
%         tpridgeTRsse = sum(trRes.^2);
%         tpridgeTRsst = sum(trY.^2);
%         tpridgeTRR2 = 1 - tpridgeTRsse/tpridgeTRsst;
%         tpridgeLen = length(ridge_selectid);
%         
%         fprintf('[ridge] train R2:%.3f\n', tpridgeTRR2);
%         
%         
%         teRes = teY - teX(:,ridge_selectid)*beta;        
%         tpridgeTEsse = sum(teRes.^2);
%         tpridgeTEsst = sum(teY.^2);
%         tpridgeTER2 = 1 - tpridgeTEsse/tpridgeTEsst;
%         fprintf('[ridge] test R2:%.3f\n', tpridgeTER2);
%     else
%         fprintf('[ridge] No var found.\n');
%     end
% end

for uu=1:MAXITER
    if mod(uu,20) == 0
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
    
    data = trX;
    response = trY;
    lambda=0.01;
    [sig_feature_l]=lasso(data, response, 'Lambda',lambda, 'Intercept',false );
    lasso_selectid = find(sig_feature_l~=0);
    
    tplassoTRsse= 0;
    tplassoTRsst = 0;
    tplassoTRR2 = 0;
    tplassoLen = 0;
    tplassoTEsse = 0;
    tplassoTEsst = 0;
    tplassoTER2 = 0;
    
    if length(lasso_selectid) > 0
        trXtruc = trX(:,lasso_selectid);
        beta=sig_feature_l(lasso_selectid);
        trRes = trY - trXtruc*beta;
        tplassoTRsse = sum(trRes.^2);
        tplassoTRsst = sum(trY.^2);
        tplassoTRR2 = 1 - tplassoTRsse/tplassoTRsst;
        tplassoLen = length(lasso_selectid);
        
        fprintf('[lasso] train R2:%.3f\n', tplassoTRR2);
        
        
        teRes = teY - teX(:,lasso_selectid)*beta;        
        tplassoTEsse = sum(teRes.^2);
        tplassoTEsst = sum(teY.^2);
        tplassoTER2 = 1 - tplassoTEsse/tplassoTEsst;
        fprintf('[lasso] test R2:%.3f\n', tplassoTER2);
    else
        fprintf('[lasso] No var found.\n');
    end
    lassoTRsse(uu,1) = tplassoTRsse;
    lassoTRsst(uu,1) = tplassoTRsst;
    lassoTRR2(uu,1) = tplassoTRR2;
    lassoLen(uu,1) = tplassoLen;
    
    lassoTEsse(uu,1) = tplassoTEsse;
    lassoTEsst(uu,1) = tplassoTEsst;
    lassoTER2(uu,1) = tplassoTER2;
end

load NIRresult01.mat
load NIRresult.mat
load NIRresult10.mat

%Q2
%mean( ridgeTER2(ridgeLen~=0) )
mean( lassoTER2(lassoLen~=0) )


%mse
%mean( ridgeTEsse(ridgeLen~=0) )
mean( lassoTEsse(lassoLen~=0) )

%how many variables found
%mean(ridgeLen(ridgeLen~=0))
mean(lassoLen(lassoLen~=0))

%how many times variables are found
%sum(ridgeLen~=0)/length(ridgeLen)
sum(lassoLen~=0)/length(lassoLen)
%% to look for lasso outcomes only
MAXITER=10
RESULTS1 = [];
RESULTS2 = [];
RESULTS3 = [];
RESULTS4 = [];
RESULTS5 = [];
for cc=1:MAXITER
    getdata = generateData();
    lambda=0.01;
    [sig_feature_l]=lasso(getdata.data, getdata.label_data, 'Lambda',lambda, 'Intercept',false);
    lasso_selectid = find(sig_feature_l~=0);

    t1 = length(lasso_selectid);
    t2 = sum( lasso_selectid < 31 );
    t3 = sum( (lasso_selectid >= 31).*(lasso_selectid < 121) );
    t4 = sum( (lasso_selectid>=121).*(lasso_selectid < 390) );
    t5 = sum( (lasso_selectid >= 391) );

    RESULTS1(cc) = t1;
    RESULTS2(cc) = t2;
    RESULTS3(cc) = t3;
    RESULTS4(cc) = t4;
    RESULTS5(cc) = t5;
    fprintf('done with OPLSoption %d, ncomp %d, cc %d\n', OPLSoption, ncomp, cc);
    fprintf('t1:%d, t2:%d, t3:%d, t4:%d, t5:%d\n', t1, t2, t3, t4,t5);
end
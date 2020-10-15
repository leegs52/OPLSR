%addpath '.\generatedata' -end
function [] = testPCRselection_newSimul()
% Not advantageous to run in a function mode. So go with a cell mode

%OPLSoptionAll = [1 2 3 4 5 6];
OPLSoptionAll = [1 2 3 5 -1];
ncompAll = [2 3 4 5];
MAXITER = 900;
RESULTS1 = zeros(length(OPLSoptionAll), length(ncompAll), MAXITER);
RESULTS2 = zeros(length(OPLSoptionAll), length(ncompAll), MAXITER);
RESULTS3 = zeros(length(OPLSoptionAll), length(ncompAll), MAXITER);
RESULTS4 = zeros(length(OPLSoptionAll), length(ncompAll), MAXITER);
RESULTS5 = zeros(length(OPLSoptionAll), length(ncompAll), MAXITER);
for bb=1:length(ncompAll)
    for cc=1:MAXITER
        getdata = generateData();
        
        for aa=1:length(OPLSoptionAll)
            OPLSoption = OPLSoptionAll(aa);
            if OPLSoption == -1
                [sig_feature] = doFDRreg(getdata.data', getdata.label_data', 0.01, [], 1);
                outcomes3.select = sig_feature(:,2);
            else
                ncomp = ncompAll(bb);
                [outcomes3] = doOneTest(getdata, OPLSoption, ncomp);
            end
            
            t1 = length(outcomes3.select);
            t2 = sum( outcomes3.select < 31 );
            t3 = sum( (outcomes3.select >= 31).*(outcomes3.select < 121) );
            t4 = sum( (outcomes3.select >= 121).*(outcomes3.select < 390) );
            t5 = sum( (outcomes3.select >= 391) );
            
            RESULTS1(aa,bb, cc) = t1;
            RESULTS2(aa,bb, cc) = t2;
            RESULTS3(aa,bb, cc) = t3;
            RESULTS4(aa,bb, cc) = t4;
            RESULTS5(aa,bb, cc) = t5;
            fprintf('done with OPLSoption %d, ncomp %d, cc %d\n', OPLSoption, ncomp, cc);
            fprintf('t1:%d, t2:%d, t3:%d, t4:%d, t5:%d\n', t1, t2, t3, t4,t5);
        end
    end
end
save('SimResult-01.mat', 'RESULTS1', 'RESULTS2', 'RESULTS3', 'RESULTS4', 'RESULTS5', 'OPLSoptionAll', 'ncompAll');

%% analysis of SimResult.mat
% (1) 1:OSC-corrected PLS RegrVector
% (2) 2:all: OSC->PLS->Reg
% (3) 3:OSC-> PLS
% (4) 5:PLS
% (5) -1:FDR
for bb=1:length(ncompAll)
    fprintf('ncomp:%d\n', ncompAll(bb));
    fprintf('t1:');
    fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        mean(RESULTS1(1,bb, :)), mean(RESULTS1(2,bb, :)), mean(RESULTS1(3,bb, :)), ...
        mean(RESULTS1(4,bb, :)), mean(RESULTS1(5,bb, :)) );
    fprintf('t2:');
    fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        mean(RESULTS2(1,bb, :)), mean(RESULTS2(2,bb, :)), mean(RESULTS2(3,bb, :)), ...
        mean(RESULTS2(4,bb, :)), mean(RESULTS2(5,bb, :)) );
    fprintf('t3:');
    fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        mean(RESULTS3(1,bb, :)), mean(RESULTS3(2,bb, :)), mean(RESULTS3(3,bb, :)), ...
        mean(RESULTS3(4,bb, :)), mean(RESULTS3(5,bb, :)) );
    fprintf('t4:');
    fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        mean(RESULTS4(1,bb, :)), mean(RESULTS4(2,bb, :)), mean(RESULTS4(3,bb, :)), ...
        mean(RESULTS4(4,bb, :)), mean(RESULTS4(5,bb, :)) );
    fprintf('t5:');
    fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        mean(RESULTS5(1,bb, :)), mean(RESULTS5(2,bb, :)), mean(RESULTS5(3,bb, :)), ...
        mean(RESULTS5(4,bb, :)), mean(RESULTS5(5,bb, :)) );
end

%% to look for FDR outcomes (-1) only
for cc=1:MAXITER
    getdata = generateData();    
    for aa=5:5
        OPLSoption = -1;
        if OPLSoption == -1
            [sig_feature] = doFDRreg(getdata.data', getdata.label_data', 0.10, [], 1);
            outcomes3.select = sig_feature(:,2);
        else
            ncomp = ncompAll(bb);
            [outcomes3] = doOneTest(getdata, OPLSoption, ncomp);
        end
        
        t1 = length(outcomes3.select);
        t2 = sum( outcomes3.select < 31 );
        t3 = sum( (outcomes3.select >= 31).*(outcomes3.select < 121) );
        t4 = sum( (outcomes3.select >= 121).*(outcomes3.select < 390) );
        t5 = sum( (outcomes3.select >= 391) );
        
        RESULTS1_5(cc) = t1;
        RESULTS2_5(cc) = t2;
        RESULTS3_5(cc) = t3;
        RESULTS4_5(cc) = t4;
        RESULTS5_5(cc) = t5;
        fprintf('done with OPLSoption %d, ncomp %d, cc %d\n', OPLSoption, ncomp, cc);
        fprintf('t1:%d, t2:%d, t3:%d, t4:%d, t5:%d\n', t1, t2, t3, t4,t5);
    end
end
save('SimResult-10-FDR.mat', 'RESULTS1_5', 'RESULTS2_5', 'RESULTS3_5', 'RESULTS4_5', 'RESULTS5_5', 'OPLSoptionAll', 'ncompAll');


%% to look for OSC-corrected PLS (1) only
ncompAll = [2 3 4 5];
alphalevels =[0.01 0.05 0.1];
RESULTS1_1 = [];
RESULTS2_1 = [];
RESULTS3_1 = [];
RESULTS4_1 = [];
RESULTS5_1 = [];
for dd=1:length(alphalevels)
    for bb=1:length(ncompAll)
        for cc=1:MAXITER
            getdata = generateData();
            
            OPLSoption = 1;
            if OPLSoption == -1
                [sig_feature] = doFDRreg(getdata.data', getdata.label_data', 0.10, [], 1);
                outcomes3.select = sig_feature(:,2);
            else
                ncomp = ncompAll(bb);
                [outcomes3] = doOneTest(getdata, OPLSoption, ncomp, alphalevels(dd));
            end
            
            t1 = length(outcomes3.select);
            t2 = sum( outcomes3.select < 31 );
            t3 = sum( (outcomes3.select >= 31).*(outcomes3.select < 121) );
            t4 = sum( (outcomes3.select >= 121).*(outcomes3.select < 390) );
            t5 = sum( (outcomes3.select >= 391) );
            
            RESULTS1_1(dd, bb, cc) = t1;
            RESULTS2_1(dd, bb, cc) = t2;
            RESULTS3_1(dd, bb, cc) = t3;
            RESULTS4_1(dd, bb, cc) = t4;
            RESULTS5_1(dd, bb, cc) = t5;
            fprintf('done with OPLSoption %d, ncomp %d, cc %d\n', OPLSoption, ncomp, cc);
            fprintf('t1:%d, t2:%d, t3:%d, t4:%d, t5:%d\n', t1, t2, t3, t4,t5);
            
        end
    end
end

save('SimResult-Num1Only.mat', 'RESULTS1_1', 'RESULTS2_1', 'RESULTS3_1', 'RESULTS4_1', 'RESULTS5_1', 'OPLSoptionAll', 'ncompAll', 'alphalevels');
%% to understand SimResult-Num1Only.mat
load SimResult-Num1Only.mat
ncompAll = [2 3 4 5];
alphalevels =[0.01 0.05 0.1];
clc
for dd=1:length(alphalevels)
    for bb=1:length(ncompAll)
        fprintf('%.3f %d:(1)%.3f (2)%.3f (3)%.3f (4)%.3f (5)%.3f\n', alphalevels(dd), ncompAll(bb), ...
            mean(RESULTS1_1(dd,bb,:)), mean(RESULTS2_1(dd,bb,:)), mean(RESULTS3_1(dd,bb,:)), ...
            mean(RESULTS4_1(dd,bb,:)), mean(RESULTS5_1(dd,bb,:)) );
    end
end


%% to check out amounts of variation
getdata = generateData();
X = getdata.data;
Y = zscore(getdata.label_data);

%I changed the following (nOrthcomp and OPTncomp)
nOrthcomp = 2;
OPTncomp=3;

meanX = mean(X,1);
meanY = mean(Y,1);
X0 = bsxfun(@minus, X, meanX);
Y0 = bsxfun(@minus, Y, meanY);

[Z,W,Pv,T] = dosc(X0,Y0,nOrthcomp,1E-3);
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(Z,Y0,OPTncomp);

TotX = sum(sum(abs(X0).^2,1));
TotZ = sum(sum(abs(Z).^2,1));
TotOrth = sum(sum(abs(T*Pv').^2,1));
TotZ + TotOrth - TotX
[TotZ  TotOrth]/TotX
PCTVAR

%% to produce regression vectors

%getdata = generateData();
%OPLSoption = 3;
OPLSoption = 1;
ncomp = 2;
[outcomes3] = doOneTest(getdata, OPLSoption, ncomp);

%% to decorate the histogram of beta_i
title('');
xlabel('Estimates of \beta_{OPLSb,i}', 'FontSize', 16);
ylabel('Frequencey', 'FontSize', 16);
%set(0, 'DefaultAxesFontSize', 16);
set(gca,'FontSize',15)

%% to decorate the regression vector
title('');
xlabel('Variable number', 'FontSize', 16);
%ylabel('OPLSa Regression Vector \beta_{OPLSa}', 'FontSize', 16);
ylabel('OPLS Regression Vector $\hat{\beta}_{OPLSb,i}$','interpreter','latex');

%set(0, 'DefaultAxesFontSize', 16);
set(gca,'FontSize',15)
%ylim([-1.7e-4 1.7e-4])
xlim([1 72]);
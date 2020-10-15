%% step 0: read data and set OPLSReg parameters
clear all; clc;
curdata = 'Neujahr_c18';

switch curdata
    case 'Neujahr_c18'
        xlsfname = 'Neujahr_c18_102011.xlsx';        
        xlssheet = 'act';
        
        dataOld = xlsread(xlsfname, xlssheet, 'E2:EB7069'); %p-by-n
        label_dataOld = xlsread(xlsfname, xlssheet, 'E7072:EB7072');%1-by-n
        [ndata, ntext, variable_nameOld] = xlsread(xlsfname, xlssheet, 'A2:A7069'); %p-by-1
        [ndata, ntext, sample_nameOld]= xlsread(xlsfname, xlssheet, 'E1:EB1');
        SampleFilter = xlsread(xlsfname, xlssheet, 'E7071:EB7071'); %using all
        VarFilter = xlsread(xlsfname, xlssheet, 'EE2:EE7069');
        responseOld = xlsread(xlsfname, xlssheet, 'E7073:EB7073');%1-by-n
        strLabelAll = {'No';'Yes'};
end

id_sample = find( SampleFilter == 1);
id_var = find( VarFilter == 1);

data = dataOld(id_var, id_sample)'; %n-by-p
if exist('label_dataOld')
    label_data = label_dataOld(1, id_sample)'; %n-by-1
end
variable_name = variable_nameOld(id_var, 1); %p-by-1
sample_name = sample_nameOld(1, id_sample)'; %n-by-1
if exist('responseOld')
    response = responseOld(1,id_sample)'; %n-by-1
end

isAutoScale = true; %it is a good idea to normalize the input responses
isFillExpect = false;

if isFillExpect
    [data] = fillDataWithExpect(data);
end
dataOrig = data;
if isAutoScale
    [data] = zscore(data); %n-by-p
end

P = size(data,2);
N = size(data,1);

if exist('label_data')
    lbl_value_all = unique(label_data)';
end


%drawfig = true;
%showSampleName = true;
%drawOneSigma = true;

%graphical setting
lw = 2;
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultAxesFontName', 'Arial');
fs = 15;
msize = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OPLSReg parameters here %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
responseNormal = zscore(response); %it is a good idea to normalize the output responses
nOrthcomp = 1; %usually, the first orthogonal component works well
alphaLevel = 0.1; %the significance level to filter variables initially using individual regression p-values.
alphaLevelP = 0.05; %the significance level to filter variables using OPLSReg vectors p-values.
%ncomp = 4;
ncomp = 3; %ncomp = nOrthcomp + # of Predictive Components. Practically, ncomp more than 2, say ncomp 3, works well.
choosetype = 'bonferoni'; % use alphaLevel/# of variables to obtain a small set of important variables
%choosetype = 'nocorrection';
saveresult = true; %to save the result of significant variables

%% step 1: selection of loadings starting with OPLS
clear params;
params.nOrthcomp = nOrthcomp;
params.ncomp = ncomp;
params.data = data;
params.responseNormal = responseNormal;

[outcomes] = doPCRselectionComps(params);
%only when you would like to override ncomp by the automatic suggestion do the following.
%ncomp = outcomes.ncomp;

%% step 2: global permutation test; proceeding with permutation tests
clear params;
params.data = data;
params.responseNormal = responseNormal;
params.alphaLevel = alphaLevel;
params.alphaLevelP = alphaLevelP;
params.choosetype = choosetype;
params.nOrthcomp = 1;
params.ncomp = ncomp;


[outcomes2] = doPCRselectionInner(params);

%% step 3: summarize the result
clear params;
params.alphaLevelP = alphaLevelP;
params.tpImportPCRcoeff = outcomes2.tpImportPCRcoeff;
params.idxAll = outcomes2.idxAll;
params.pcr_vector = outcomes2.pcr_vector;
params.RegPvalues = outcomes2.RegPvalues;
params.RegCoeffs = outcomes2.RegCoeffs;
params.doSave = saveresult;

outcomes3 = doPCRselectionSummarize(params);
% %% fdr, lasso
% lambda=0.25;
% [sig_feature_l]=lasso(data, response, 'Lambda',lambda, 'Intercept',false);
% lasso_selectid = find(sig_feature_l~=0);
% length(lasso_selectid)
% 
% alphaLevel = 0.000000000001;
% qlev = alphaLevel;
% 
% chem_shift_FDR = (1:size(data,2))'; %p-by-1
% [sig_feature] = doFDRreg(data', response', qlev, chem_shift_FDR, 1);
% fdrselectid = sig_feature(:,end);
% length(fdrselectid)
%% analyze the outcome

function [outcomes3] = doOneTest(getdata, OPLSoption, ncomp, alphaLevelParam)
%getdata = generateData();
data = getdata.data;
response = getdata.label_data;

responseNormal = zscore(response);
nOrthcomp = 1;
%ncomp = 8;
%alphaLevel = 0.05;
%alphaLevel = 0.10;
%alphaLevel = 0.01;
if nargin < 4
    alphaLevel = 0.1;
else
    alphaLevel = alphaLevelParam;
end

alphaLevelP = 0.05;


%% step 1: selection of loadings starting with OPLS
% clear params;
% params.nOrthcomp = nOrthcomp;
% params.ncomp = ncomp;
% params.data = data;
% params.responseNormal = responseNormal;
% [outcomes] = doPCRselectionComps(params);
% OPTtype = outcomes.OPTtype;
% idPCAll = outcomes.idPCAll;
% ncomp = outcomes.ncomp;

%% step 2: global permutation test; proceeding with permutation tests
clear params;
params.data = data;
params.responseNormal = responseNormal;
params.alphaLevel = alphaLevel;
params.alphaLevelP = alphaLevelP;
%params.choosetype = 'bonferoni';
params.choosetype = 'nocorrection';
params.nOrthcomp = 1;
params.ncomp = ncomp;
params.OPTtype = 'OPLS';
%params.idPCAll = idPCAll;
params.OPLSoption = OPLSoption;
%params.minITERset = 30;
%params.maxITERset = 300;

[outcomes2] = doPCRselectionInner(params);

%% summarize the result
clear params;
params.alphaLevelP = alphaLevelP;
params.tpImportPCRcoeff = outcomes2.tpImportPCRcoeff;
params.idxAll = outcomes2.idxAll;
params.pcr_vector = outcomes2.pcr_vector;
params.RegPvalues = outcomes2.RegPvalues;
params.RegCoeffs = outcomes2.RegCoeffs;
params.doSave = false;

outcomes3 = doPCRselectionSummarize(params);
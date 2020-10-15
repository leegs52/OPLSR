function [outcomes] = doPCRselectionSummarize(params)


if isfield(params, 'verbose')
    verbose = params.verbose;
else
    verbose = true;
    %OPLSIncludeOrtho = 0;
end

alphaLevelP = params.alphaLevelP;

tpImportPCRcoeff = params.tpImportPCRcoeff;
idxAll = params.idxAll;
pcr_vector = params.pcr_vector;
RegPvalues = params.RegPvalues;
RegCoeffs = params.RegCoeffs;
doSave = params.doSave;

%doIndCutoff = false;
doIndCutoff = true;
if length(params.tpImportPCRcoeff) == 0
    outcomes.outcomeArray = [];
    outcomes.select = [];
    return;
end

[IDfoundPCR, lwcutoffA, upcutoffA, achvSigLevel, tpZeroAll, tpdtrankaSigLevel] = innerSummarize(pcr_vector, idxAll, tpImportPCRcoeff, alphaLevelP, doIndCutoff, verbose);

% % individual, for doIndCutoff == true
% lwcutoffAInd = zeros(length(idxAll),1);
% upcutoffAInd = zeros(length(idxAll),1);
% for ee=1:length(idxAll)
%     lwcutoff = quantile(tpImportPCRcoeff(ee).coeff,alphaLevelP/2);
%     upcutoff =quantile(tpImportPCRcoeff(ee).coeff,1-alphaLevelP/2);
%     lwcutoffAInd(ee,1) = lwcutoff;
%     upcutoffAInd(ee,1) = upcutoff;
% end
% 
% %prepare the global quantiles
% tpAllcoeff = [];
% for kk=1:length(tpImportPCRcoeff)
%     tpAllcoeff = [tpAllcoeff; tpImportPCRcoeff(kk).coeff];
% end
% 
% if verbose
%     figure(35); clf; hold on;
%     hist(tpAllcoeff, min(200, length(tpAllcoeff)/100));
%     title('histogram of pcr coefficients when variables do not relate to Y');
% end
% 
% % [hdecision, hpval] = jbtest(tpAllcoeff);
% % if hdecision == 1
% %     title('The density is not from a normal distribution');
% % else
% %     title('The density is from a normal distribution')
% % end
% 
% lwcutoffAll = quantile(tpAllcoeff,alphaLevelP/2);
% upcutoffAll =quantile(tpAllcoeff,1-alphaLevelP/2);
% lwcutoffAAll = zeros(length(idxAll),1);
% upcutoffAAll = zeros(length(idxAll),1);
% lwcutoffAAll(1:length(idxAll),1) = lwcutoffAll;
% upcutoffAAll(1:length(idxAll),1) = upcutoffAll;
% 
% 
% %fprintf('[global]\t %.3e\t %.3e \n',lwcutoffAll, upcutoffAll);
% %fprintf('[ind_mean]\t %.3e\t %.3e \n',mean(lwcutoffAInd), mean(upcutoffAInd));
% 
% IDfoundPCR = zeros(length(idxAll),1);
% if doIndCutoff
%     lwcutoffA = lwcutoffAInd;
%     upcutoffA = upcutoffAInd;
% else
%     lwcutoffA = lwcutoffAAll;
%     upcutoffA = upcutoffAAll;
% end
% 
% achvSigLevel =  zeros(length(idxAll),1);
% for ee=1:length(idxAll)
%     tplwcutoff = quantile(tpImportPCRcoeff(ee).coeff,alphaLevelP/2);
%     tpupcutoff =quantile(tpImportPCRcoeff(ee).coeff,1-alphaLevelP/2);    
%     
%     %achvSigLevel(ee,1) = 2*sum(tpImportPCRcoeff(ee).coeff > abs(pcr_vector(idxAll(ee))))/length(tpImportPCRcoeff(ee).coeff);
%     achvSigLevel(ee,1) = (sum(tpImportPCRcoeff(ee).coeff > abs(pcr_vector(idxAll(ee)))) + ...
%         sum(tpImportPCRcoeff(ee).coeff < -abs(pcr_vector(idxAll(ee)))) )...
%         /length(tpImportPCRcoeff(ee).coeff);
%     
%     %fprintf('Now %d. The # of: %d. The achivSig: %.4f \n', idxAll(ee), length(tpImportPCRcoeff(ee).coeff), achvSigLevel(ee,1));
%     if achvSigLevel(ee,1) < alphaLevelP
%         %if ( abs(pcr_vector(idxAll(ee)))-tplwcutoff) * ( abs(pcr_vector(idxAll(ee))) - tpupcutoff ) < 0
%         %    fprintf('Strange!\n');
%         %end
%         if verbose
%             %fprintf('Found %d at %d The # of permuted coefficients: %d. The achivSig: %.4f \n', idxAll(ee), ee, length(tpImportPCRcoeff(ee).coeff), achvSigLevel(ee,1));
%         end
%         IDfoundPCR(ee)=1;
%     end
% end
% 
% if verbose
%     fprintf('Overall # variables found:%d\n', sum(IDfoundPCR));
% end

% figure(34); clf; hold on;
% plot(pcr_vector, 'b.');
% %use bonferoni for this time. It should be set by the policy.
% %Look at idxAll also.
% %tpIDreg = idxAll;
% %plot(idxAll, pcr_vector(idxAll), 'o', 'Color', 'g');
% plot([1 length(pcr_vector)], [0 0], ':', 'Color', 'k');
% %plot(idxAll, RegCoeffs(idxAll)/max(RegCoeffs(idxAll))*max(pcr_vector(idxAll)), '*', 'LineWidth', 2, 'Color', 'c');
% 
% %[AX,H1,H2] = plotyy(idxAll, pcr_vector(idxAll),idxAll, RegPvalues(idxAll),'plot');
% %set(H1,'LineStyle','og')
% %set(H2,'LineStyle','sg')
% for kk=1:length(idxAll)
%    plot(idxAll(kk), lwcutoffA(kk) , 'rv', 'linewidth', 1, 'MarkerFaceColor', 'r');
%    plot(idxAll(kk), upcutoffA(kk) , 'r^', 'linewidth', 1, 'MarkerFaceColor', 'r');
%    plot([idxAll(kk) idxAll(kk)] , [lwcutoffA(kk) upcutoffA(kk)] , ':r', 'linewidth', 1);
% end
% %plot(idxAll, pcr_vector(idxAll), 'o', 'Color', 'g');
% plot(idxAll, pcr_vector(idxAll), 'o', 'Color', 'g', 'linewidth', 3);
% xlim([1 length(pcr_vector)]);
% ylabel('OPLS Regression Vector');
% xlabel('Variables');
% %plot([idxAll idxAll], [lwcutoffA lwcutoffA], '-r');
% %plot([idxAll idxAll], [upcutoffA upcutoffA], '-r');
% 
% 
% %gather information on orders. Is it necessary? Maybe we don't need it.
% tpdt = pcr_vector(idxAll);
% [tpdtsort tpidxs] = sort(achvSigLevel , 'ascend');
% tpdtrankaSigLevel = zeros(length(idxAll),1);
% tpdtrankaSigLevel(tpidxs) = (1:length(idxAll))';
% 
% tpZeroAll = zeros(length(idxAll),1);
% tpZeroAll(find(IDfoundPCR==1)) = 1; %the variables that were significant by permutation
% 
% %plot(idxAll(tpZeroAll==1), pcr_vector(idxAll(tpZeroAll==1)), 'd', 'Color', 'k', 'linewidth', 3);
% plot(idxAll(tpZeroAll==1), pcr_vector(idxAll(tpZeroAll==1)), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'linewidth', 3);




%1:variable-ID
%2:individual pvalue for y ~ individual x
%3:individual regression coefficient for y ~ individual x
%4:1 if significant by permuation, 0 if not.
%5:acheived signifiance level from permutation
%6:the rank of acheived signifiance levels
%7:coefficient PC regression vector

outcomeArray = [idxAll RegPvalues(idxAll) RegCoeffs(idxAll) tpZeroAll achvSigLevel tpdtrankaSigLevel pcr_vector(idxAll)];
outcomes.outcomeArray = outcomeArray;
outcomes.select = idxAll(tpZeroAll==1);
if doSave
    csvwrite('tpPCRoutput.csv',outcomeArray);
end

%additional weight vector analysis
%[IDfoundPCRw, lwcutoffAw, upcutoffAw, achvSigLevelw, tpZeroAllw, tpdtrankaSigLevelw] = innerSummarize(params.OrthWVector, idxAll, params.tpOrthWVector, alphaLevelP, doIndCutoff, verbose);




function [IDfoundPCR, lwcutoffA, upcutoffA, achvSigLevel, tpZeroAll, tpdtrankaSigLevel] = innerSummarize(pcr_vector, idxAll, tpImportPCRcoeff, alphaLevelP, doIndCutoff, verbose)

% individual, for doIndCutoff == true
lwcutoffAInd = zeros(length(idxAll),1);
upcutoffAInd = zeros(length(idxAll),1);
for ee=1:length(idxAll)
    lwcutoff = quantile(tpImportPCRcoeff(ee).coeff,alphaLevelP/2);
    upcutoff =quantile(tpImportPCRcoeff(ee).coeff,1-alphaLevelP/2);
    lwcutoffAInd(ee,1) = lwcutoff;
    upcutoffAInd(ee,1) = upcutoff;
end

%additional distribution
%eee=find(upcutoffAInd-lwcutoffAInd==max(upcutoffAInd-lwcutoffAInd))
%hist(tpImportPCRcoeff(eee).coeff,100)

%prepare the global quantiles
tpAllcoeff = [];
for kk=1:length(tpImportPCRcoeff)
    tpAllcoeff = [tpAllcoeff; tpImportPCRcoeff(kk).coeff];
end

if verbose
    figure(35); clf; hold on;
    hist(tpAllcoeff, min(200, length(tpAllcoeff)/100));
    title('histogram of pcr coefficients when variables do not relate to Y');
end

% [hdecision, hpval] = jbtest(tpAllcoeff);
% if hdecision == 1
%     title('The density is not from a normal distribution');
% else
%     title('The density is from a normal distribution')
% end

lwcutoffAll = quantile(tpAllcoeff,alphaLevelP/2);
upcutoffAll =quantile(tpAllcoeff,1-alphaLevelP/2);
lwcutoffAAll = zeros(length(idxAll),1);
upcutoffAAll = zeros(length(idxAll),1);
lwcutoffAAll(1:length(idxAll),1) = lwcutoffAll;
upcutoffAAll(1:length(idxAll),1) = upcutoffAll;


%fprintf('[global]\t %.3e\t %.3e \n',lwcutoffAll, upcutoffAll);
%fprintf('[ind_mean]\t %.3e\t %.3e \n',mean(lwcutoffAInd), mean(upcutoffAInd));

IDfoundPCR = zeros(length(idxAll),1);
if doIndCutoff
    lwcutoffA = lwcutoffAInd;
    upcutoffA = upcutoffAInd;
else
    lwcutoffA = lwcutoffAAll;
    upcutoffA = upcutoffAAll;
end

achvSigLevel =  zeros(length(idxAll),1);
for ee=1:length(idxAll)
    tplwcutoff = quantile(tpImportPCRcoeff(ee).coeff,alphaLevelP/2);
    tpupcutoff =quantile(tpImportPCRcoeff(ee).coeff,1-alphaLevelP/2);    
    
    %achvSigLevel(ee,1) = 2*sum(tpImportPCRcoeff(ee).coeff > abs(pcr_vector(idxAll(ee))))/length(tpImportPCRcoeff(ee).coeff);
    %Change the application of 
    achvSigLevel(ee,1) = (sum(tpImportPCRcoeff(ee).coeff > abs(pcr_vector(idxAll(ee)))) + ...
        sum(tpImportPCRcoeff(ee).coeff < -abs(pcr_vector(idxAll(ee)))) )...
        /length(tpImportPCRcoeff(ee).coeff);
    
    %fprintf('Now %d. The # of: %d. The achivSig: %.4f \n', idxAll(ee), length(tpImportPCRcoeff(ee).coeff), achvSigLevel(ee,1));
    %Change the application of critical values on April 29, 2012
    %if achvSigLevel(ee,1) < alphaLevelP
    %if ( abs(pcr_vector(idxAll(ee)))-tplwcutoff) * ( abs(pcr_vector(idxAll(ee))) - tpupcutoff ) < 0
    if ( abs(pcr_vector(idxAll(ee)))-tplwcutoff) * ( abs(pcr_vector(idxAll(ee))) - tpupcutoff ) > 0
        %    fprintf('Strange!\n');
        %end
        if verbose
            %fprintf('Found %d at %d The # of permuted coefficients: %d. The achivSig: %.4f \n', idxAll(ee), ee, length(tpImportPCRcoeff(ee).coeff), achvSigLevel(ee,1));
        end
        IDfoundPCR(ee)=1;
    end
end

if verbose
    fprintf('Overall # variables found:%d\n', sum(IDfoundPCR));
end

figure(34); clf; hold on;
plot(pcr_vector, 'b.');
%use bonferoni for this time. It should be set by the policy.
%Look at idxAll also.
%tpIDreg = idxAll;
%plot(idxAll, pcr_vector(idxAll), 'o', 'Color', 'g');
plot([1 length(pcr_vector)], [0 0], ':', 'Color', 'k');
%plot(idxAll, RegCoeffs(idxAll)/max(RegCoeffs(idxAll))*max(pcr_vector(idxAll)), '*', 'LineWidth', 2, 'Color', 'c');

%[AX,H1,H2] = plotyy(idxAll, pcr_vector(idxAll),idxAll, RegPvalues(idxAll),'plot');
%set(H1,'LineStyle','og')
%set(H2,'LineStyle','sg')
for kk=1:length(idxAll)
   plot(idxAll(kk), lwcutoffA(kk) , 'rv', 'linewidth', 1, 'MarkerFaceColor', 'r','MarkerSize', 5);
   plot(idxAll(kk), upcutoffA(kk) , 'r^', 'linewidth', 1, 'MarkerFaceColor', 'r','MarkerSize', 5);
   plot([idxAll(kk) idxAll(kk)] , [lwcutoffA(kk) upcutoffA(kk)] , ':r', 'linewidth', 1);
end
%plot(idxAll, pcr_vector(idxAll), 'o', 'Color', 'g');
plot(idxAll, pcr_vector(idxAll), 's', 'Color', 'g', 'linewidth', 2);
xlim([1 length(pcr_vector)]);
ylabel('OPLS Regression Vector $\hat{\beta}_{OPLSb,i}$','interpreter','latex');
xlabel('Variable number');
%plot([idxAll idxAll], [lwcutoffA lwcutoffA], '-r');
%plot([idxAll idxAll], [upcutoffA upcutoffA], '-r');


%gather information on orders. Is it necessary? Maybe we don't need it.
tpdt = pcr_vector(idxAll);
[tpdtsort tpidxs] = sort(achvSigLevel , 'ascend');
tpdtrankaSigLevel = zeros(length(idxAll),1);
tpdtrankaSigLevel(tpidxs) = (1:length(idxAll))';

tpZeroAll = zeros(length(idxAll),1);
tpZeroAll(find(IDfoundPCR==1)) = 1; %the variables that were significant by permutation

%plot(idxAll(tpZeroAll==1), pcr_vector(idxAll(tpZeroAll==1)), 'd', 'Color', 'k', 'linewidth', 3);
plot(idxAll(tpZeroAll==1), pcr_vector(idxAll(tpZeroAll==1)), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'linewidth', 2);
%sigft=idxAll(tpZeroAll==1)




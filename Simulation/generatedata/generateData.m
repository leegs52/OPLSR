function [gendata]=generateData()
%clear all; clc;
doVerbose = false;
N = 60;
%P = 400;
P = 1000;
tpRandnum = rand(N,1);
label_data = (tpRandnum > 1/3) + (tpRandnum > 2/3) + 1;

data = zeros(N, P);

idx_1 = find( label_data == 1);
idx_2 = find( label_data == 2);
idx_3 = find( label_data == 3);

for kk=1:P
    if mod(kk,300) == 0
        if doVerbose
            fprintf('Now kk %d\n',kk);
        end
    end
    if kk >= 1 && kk <= 4
        data(idx_1, kk) = 5 + rand( length(idx_1), 1);
        data(idx_2, kk) = 3 + rand( length(idx_2), 1);
        data(idx_3, kk) = 1 + rand( length(idx_3), 1);
    elseif kk >= 5 && kk <= 8
        data(idx_1, kk) = 1 + rand( length(idx_1), 1);
        data(idx_2, kk) = 3 + rand( length(idx_2), 1);
        data(idx_3, kk) = 5 + rand( length(idx_3), 1);
    elseif kk== 9 || kk == 11 || kk == 13 || kk == 15 || kk == 29
        [r r2 r3] = gen2Dsample3(length(idx_1), length(idx_2), length(idx_3));
        data(idx_1, kk) = r(:,1);
        data(idx_1, kk+1) = r(:,2);
        data(idx_2, kk) = r2(:,1);
        data(idx_2, kk+1) = r2(:,2);
        data(idx_3, kk) = r3(:,1);
        data(idx_3, kk+1) = r3(:,2);
    elseif kk == 17 || kk == 20 || kk == 23 || kk == 26
        [r1 r2 r3] = gen3Dsample3(length(idx_1), length(idx_2), length(idx_3));
        data(idx_1, kk) = r1(:,1);
        data(idx_1, kk+1) = r1(:,2);
        data(idx_1, kk+2) = r1(:,3);
        data(idx_2, kk) = r2(:,1);
        data(idx_2, kk+1) = r2(:,2);
        data(idx_2, kk+2) = r2(:,3);
        data(idx_3, kk) = r3(:,1);
        data(idx_3, kk+1) = r3(:,2);
        data(idx_3, kk+2) = r3(:,3);
    elseif ( (kk >= 31) && (kk < 121))
        if mod(kk-31,3) == 0
            for ss=1:size(data,1)
                data(ss, kk:(kk+2))=gennetvar(data( ss, floor((kk-30)/3)+1), 10, 3);
            end
        end
    elseif ( (kk >= 121) && (kk < 391))
        if mod(kk-121,3) == 0
            for ss=1:size(data,1)
                data(ss, kk:(kk+2))=gennetvar(data( ss, floor((kk-120)/3)+30), 10, 3);
            end
        end
    elseif ( (kk >= 391))
        data(idx_1, kk) = 2 + rand( length(idx_1), 1);
        data(idx_2, kk) = 2 + rand( length(idx_2), 1);
        data(idx_3, kk) = 2 + rand( length(idx_3), 1);
        [tb,tbint,tr,trint,tstats] = regress(label_data,[ones(size(data,1),1)  data(:,kk)]);
        while tstats(3) < 0.1
            data(idx_1, kk) = 2 + rand( length(idx_1), 1);
            data(idx_2, kk) = 2 + rand( length(idx_2), 1);
            [tb,tbint,tr,trint,tstats] = regress(label_data,[ones(size(data,1),1)  data(:,kk)]);
        end
    end
end

gendata.data = data;
gendata.idx_1 = idx_1;
gendata.idx_2 = idx_2;
gendata.idx_3 = idx_3;
gendata.label_data = label_data;

doShowThree = true;
if doShowThree
    figure(54); clf; hold on;
    kk=17; mksize = 13; lw=2;
    plot3(data(idx_1, kk), data(idx_1, kk+1), data(idx_1, kk+2), 'bo', 'LineWidth', lw, 'MarkerSize', mksize);
    plot3(data(idx_2, kk), data(idx_2, kk+1), data(idx_2, kk+2), 'rx', 'LineWidth', lw, 'MarkerSize', mksize);
    plot3(data(idx_3, kk), data(idx_3, kk+1), data(idx_3, kk+2), 'gs', 'LineWidth', lw, 'MarkerSize', mksize);
    legend('Label 1', 'Label 2', 'Label 3');
    grid on;
    xlabel('X_1');ylabel('X_2');zlabel('X_3');
    view([46 32]);
    %saveas(gcf, 'illust_X1X2X3.fig');
    %testt(data(idx_1, kk), data(idx_2, kk))
    %testt(data(idx_1, kk+1), data(idx_2, kk+1))
    %testt(data(idx_1, kk+2), data(idx_2, kk+2))
end

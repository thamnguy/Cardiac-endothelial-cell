uniqueSampleID = {'Fetal', 'CTL-P1', 'CTL-P28', 'CTL-P56', '', 'AR(1)-P28', 'AR(1)-P56','',  ...
    'MI(28)-P30', 'AR(1)+MI(28)-P30', '', 'MI(28)-P35', 'AR(1)+MI(28)-P35', '',  ...
    'MI(28)-P42', 'AR(1)+MI(28)-P42', '', 'MI(28)-P56', 'AR(1)+MI(28)-P56'};
tickLbl = {'Fetal', 'CTL-P1', 'CTL-P28', 'CTL-P56', '', 'AR_P_1P28', 'AR_P_1P56','',  ...
    'MI_P_2_8P30', 'AR_P_1MI_P_2_8P30', '', 'MI_P_2_8P35', 'AR_P_1MI_P_2_8P35', '',  ...
    'MI_P_2_8P42', 'AR_P_1MI_P_2_8P42', '', 'MI_P_2_8P56', 'AR_P_1MI_P_2_8P56'};
countCell = zeros(max(idx), length(uniqueSampleID));
for i = 1 : max(idx)
    for j = 1 : length(uniqueSampleID)
        countCell(i, j) = length( find( idx==i & ismember(sampleID, uniqueSampleID(j)) ) );
    end
end

ratioMat = zeros(size(countCell));
for i = 1 :  length(uniqueSampleID)
    ratioMat(:, i) = countCell(:, i) / sum(countCell(:, i)) *100;
end

bar(ratioMat','stacked');
legendLbl = {'VEC1', 'VEC2', 'VEC3', 'LEC1', 'LEC2'};
legend(legendLbl, 'Location','eastoutside'); 
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 16)
ylabel('%EC');
xticks(1:1:length(uniqueSampleID)); xticklabels(tickLbl);
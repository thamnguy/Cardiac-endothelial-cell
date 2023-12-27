uniqueSampleID = {'Fetal', 'CTL-P1', 'CTL-P28', 'CTL-P56', '', 'AR(1)-P28', 'AR(1)-P56','',  ...
    'MI(28)-P30', 'AR(1)+MI(28)-P30', '', 'MI(28)-P35', 'AR(1)+MI(28)-P35', '',  ...
    'MI(28)-P42', 'AR(1)+MI(28)-P42', '', 'MI(28)-P56', 'AR(1)+MI(28)-P56'};
tickLbl = {'Fetal', 'CTL-P1', 'CTL-P28', 'CTL-P56', '', 'AR_P_1P28', 'AR_P_1P56','',  ...
    'MI_P_2_8P30', 'AR_P_1MI_P_2_8P30', '', 'MI_P_2_8P35', 'AR_P_1MI_P_2_8P35', '',  ...
    'MI_P_2_8P42', 'AR_P_1MI_P_2_8P42', '', 'MI_P_2_8P56', 'AR_P_1MI_P_2_8P56'};

sampleGroup = { {'FH1', 'FH2', 'FH3'}, {'8026-P1', '8094-P1', '8095-P1'}, ...
    {'8046-BZ', '8231-BZ', '8232-BZ'}, {'8052-AZ', '8233-BZ', '8234-BZ'}, {}, ...
    {'8030-CZ', '8014-BZ', '8015-BZ', '8015-AZ'}, {'8052-AZ', '8233-BZ', '8234-BZ'}, {}, ...
    {'8026-BZ', '8178-BZ'}, {'8064-AZ', '8064-CZ', '8195-BZ', '8213-BZ', '8216-BZ'}, {}, ...
    {'8179-BZ', '8205-BZ'}, {'8065-AZ', '8065-CZ', '8095-AZ', '8095-BZ', '8194-BZ'}, {}, ...
    {'8210-BZ', '8211-BZ'}, {'8094-AZ', '8215-BZ', '8217-BZ'}, {}, ...
    {'8085-BZ', '8086-BZ'}, {'7995-BZ', '8060-AZ', '8060-IZ'}, ...
    };

%% VEC1
figure, hold on
countCell = zeros(1, length(uniqueSampleID));
for j = 1 : length(uniqueSampleID)
    countCell(1, j) = length( find( idx==1 & ismember(sampleID, uniqueSampleID(j)) == 1 ) ) / ...
        length( find(ismember(sampleID, uniqueSampleID(j)) == 1) ) * 100;
end
bar(1:length(uniqueSampleID), countCell, 'FaceColor', 'red'); 

for i = 1 : length(uniqueSampleID)
    theGroup = sampleGroup{i};
    for j = 1 : length(theGroup)
        percenTage = length( find( idx==1 & ismember(originalSampleID, theGroup(j)) == 1 ) ) / ...
        length( find(ismember(originalSampleID, theGroup(j)) == 1) ) * 100;
         plot(i, percenTage, 'Marker', '.', 'Color', 'k', 'MarkerSize', 15, 'LineStyle', 'none')
    end
end
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 16)
ylabel('%EC'); title('VEC1')
xticks(1:1:length(uniqueSampleID)); xticklabels(tickLbl);
hold off

%% VEC3
figure, hold on
countCell = zeros(1, length(uniqueSampleID));
for j = 1 : length(uniqueSampleID)
    countCell(1, j) = length( find( idx==3 & ismember(sampleID, uniqueSampleID(j)) == 1 ) ) / ...
        length( find(ismember(sampleID, uniqueSampleID(j)) == 1) ) * 100;
end
bar(1:length(uniqueSampleID), countCell, 'FaceColor', 'green'); 

for i = 1 : length(uniqueSampleID)
    theGroup = sampleGroup{i};
    for j = 1 : length(theGroup)
        percenTage = length( find( idx==3 & ismember(originalSampleID, theGroup(j)) == 1 ) ) / ...
        length( find(ismember(originalSampleID, theGroup(j)) == 1) ) * 100;
         plot(i, percenTage, 'Marker', '.', 'Color', 'k', 'MarkerSize', 15, 'LineStyle', 'none')
    end
end
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 16)
ylabel('%EC'); title('VEC3')
xticks(1:1:length(uniqueSampleID)); xticklabels(tickLbl);
hold off

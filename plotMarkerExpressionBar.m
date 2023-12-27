function [] = plotMarkerExpressionBar(marker, gene, expressMat, groupID)

if ~exist('gene','var')
    error('\r\nThe gene list (name: gene) has not been loaded or aciddentally deleted\r\n')
    return;
end

[~, index] = ismember(marker, gene);
if index < 1
    error('\r\nThe marker gene does not exist in the gene list\r\n')
    return;
end

if ~exist('expressMat','var')
    error('\r\nThe expression matrix (name: expressMat) has not been loaded or aciddentally deleted\r\n')
    return;
end

markerExpression = full( expressMat(index, :) );
uniqueSampleID = unique(groupID);

meanExp = zeros(size(uniqueSampleID));
stdExp = zeros(size(uniqueSampleID));
for i = 1 : length(uniqueSampleID)
    meanExp(i) = mean( markerExpression( ismember(groupID, uniqueSampleID(i)) ) );
    stdExp(i) = std( markerExpression( ismember(groupID, uniqueSampleID(i)) ) ) ...
        / sqrt( length( find( ismember(groupID, uniqueSampleID(i)) == 1) ) );
end
figure, bar(1:length(uniqueSampleID), meanExp);
hold on

er = errorbar(1:length(uniqueSampleID), meanExp, stdExp/1, stdExp/1);
er.Color = [0 0 0];
er.LineStyle = 'none';

set(gca,'FontSize',16)
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off')
ylabel( 'average expression' );
xticks(1:1:length(uniqueSampleID));xticklabels(uniqueSampleID);

title(marker);
hold off


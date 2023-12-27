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
figure, violinplot(markerExpression, groupID)

set(gca,'FontSize',16)
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off')
ylabel( ' expression' );

title(marker);
hold off


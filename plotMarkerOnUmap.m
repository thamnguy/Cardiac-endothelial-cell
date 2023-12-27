function [] = plotMarkerOnUmap(marker, threshold, gene, expressMat, coordinate)
% This function plot the cells expressing the specific genes in the 2D
% (localization plot). The expression level is reflected by heat-map color.
% marker: gene name to be plot
% gene: all gene names in the data
% expressMat: the expression matrix. Each row corresponds to a gene. The
%             number of (all) genes is the number of row. Each cell
%             corresponds to a cell.
% coordinate: the 2D coordinate. Each row correspond to a cells. The first
%             column is the X-coordinate of the cell on the 2D
%             visualization. The second column is the Y-coordinate.
% threshold: the minimum expression level of the cells that can be
%             reflected by heatmap color. Cells that have expression <=
%             threshold would be represented by light gray.

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

if ~exist('coordinate','var')
    error('\r\nThe umap (name: coordinate) has not been loaded or aciddentally deleted\r\n')
    return;
end
if length(markerExpression) ~= length(coordinate(:, 1))
    error('\r\nThe umap and expression data are not matched\r\n')
    return;
end

index1 = find(markerExpression <= threshold );
figure, scatter(coordinate(index1, 1), coordinate(index1, 2), 15, [192, 192, 192]/255, '.');
hold on
index1 = find(markerExpression > threshold );
scatter(coordinate(index1, 1), coordinate(index1, 2), 15, markerExpression(index1), '.');
colorbar;
colormap('jet');
title([marker, ' > ', num2str(threshold)])
xlabel('umap 1'); ylabel('umap 2');
% title([ marker,  ' > ', num2str(threshold)]);
set(gca,'FontSize',16)
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off')
hold off


logScaleExpress = log2(expressMat + 1);

allECMarkerJoin = { 'COL1A1', 'COL1A2', 'FN1', 'PDGFRA',...
    'ACVR1', 'LTBP1', 'NEO1', 'BMPR1A', ...
    'EXT1', 'HIVEP1', 'SMAD6', 'BMP6' };
ylabel = {'COL1A1', 'COL1A2', 'FN1', 'PDGFRA',...
    'ACVR1', 'LTBP1', 'NEO1', 'BMPR1A', ...
    'EXT1', 'HIVEP1', 'SMAD6', 'BMP6'};


[~, geneIndex] = ismember( allECMarkerJoin, gene );
plotLogScaleExpress = logScaleExpress(geneIndex, :);

meanValList = zeros(size(plotLogScaleExpress, 1), 1);
stdValList = zeros(size(plotLogScaleExpress, 1), 1);
for i = 1 : size(plotLogScaleExpress, 1)
    meanValList(i) = mean(plotLogScaleExpress(i, :));
    stdValList(i) = std(plotLogScaleExpress(i, :));
end

plotMeanExp = zeros(length(geneIndex), max(idx));
for i = 1 : length(geneIndex)
    for j = 1 : max(idx)
        cellIndex = find(idx==j);
        plotMeanExp(i, j) = ( mean(plotLogScaleExpress(i, cellIndex)) - meanValList(i) ) / stdValList(i);
    end
end

figure, heatmap(plotMeanExp); colormap('jet');
Ax = gca;
Ax.XDisplayLabels = nan( size(Ax.XDisplayData) );
%Ax.YDisplayLabels = nan( size(Ax.YDisplayData) );
Ax.YDisplayLabels = ylabel;
Ax.XDisplayLabels = {'VEC1', 'VEC2', 'VEC3', 'LEC1', 'LEC2'};
Ax.ColorLimits = [-0.5 1.5];
Ax.CellLabelColor = 'none';
set(gca, 'FontName', 'Arial', 'FontSize', 16)

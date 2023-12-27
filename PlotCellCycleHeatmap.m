logScaleExpress = log2(expressMat + 1);

%% fold - change heatmap (Z-score) of cell cycle
allECMarkerJoin = { 'AURKB', 'ENSSSCG00000026302' 'INCENP', 'CDCA8', 'BIRC5', 'POLA1', ...
    'POLE', 'DNA2', 'ENSSSCG00000026055', 'DDB2', 'FANCI', 'FANCM', 'TOP2A', 'HMGB2', 'NCAPG', ...
    'E2F8', 'FOXM1', 'E2F7', 'MYBL1', 'E2F1', 'MXD3', 'MYBL2', 'PBX4', 'DNMT1', 'ARNTL2', ...
    'LEF1', 'EBF4', 'E2F3', 'TFDP1', 'LIN54', 'ETV4', 'PRDM5', 'MBD3', 'MAZ', 'MAFG', 'TERF1', 'PBX3'};
ylabel = {'AURKB', 'MKI67' 'INCENP', 'CDCA8', 'BIRC5', 'POLA1', ...
    'POLE', 'DNA2', 'PRIM1', 'DDB2', 'FANCI', 'FANCM', 'TOP2A', 'HMGB2', 'NCAPG', ...
    'E2F8', 'FOXM1', 'E2F7', 'MYBL1', 'E2F1', 'MXD3', 'MYBL2', 'PBX4', 'DNMT1', 'ARNTL2', ...
    'LEF1', 'EBF4', 'E2F3', 'TFDP1', 'LIN54', 'ETV4', 'PRDM5', 'MBD3', 'MAZ', 'MAFG', 'TERF1', 'PBX3'};


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

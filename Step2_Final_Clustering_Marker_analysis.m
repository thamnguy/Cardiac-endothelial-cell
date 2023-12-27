%% Manuscript: Nguyen et al. Single-cell RNA sequencing analysis identifies 
% one subpopulation of endothelial cells that proliferates and another that 
% undergoes the endothelial-mesenchymal transition in regenerating pig
% hearts. Frontiers in Bioengineering and Biotechnology Tissue Engineering 
% and Regenerative Medicine. 2023. DOI: 10.3389/fbioe.2023.1257669. This
% manuscript examine endothelial cells and angiogenesis in pig hearts that
% recover from myocardial infarction induced 4 weeks after birth.

% This is Step 2 of the endothelial cells single nuclei RNA sequencing (snRNAseq)
% analyusis presented in the manuscript above. After this step, endothelial
% cells are grouped into five supopulations; including one subpopulation of
% proliferating endothelial cell (indication of angiogenesis) and one
% subpopulation of endothelial cell expressing endo-mesenchymal cluster,
% which only appears in recovered hearts.

% Results from 'Step 1' are already precomputed, so Step 2 can be executed
% without reruning Step 1. If rerunning step 1, be careful in cluster
% number assignment (see below).

%% clear the workspace (if needed)
clear all
close all
clc

%% load the endothelial cell pre-cluster (step 1) data
% load the endothelial cell snRNAseq data. The data is already normalized in 'step 0' 
load EC_expressMat.mat 

% load the entire gene list in the data
load allGene.mat;

% load the UMAP and pre-cluster result from step 1
load allEC_GOCycle_Umap.mat % umap 2D coordinate
load allEC_GOCycle_UmapClusterID.mat % umap precluster

% load the pig group for each cell
load cell_pigGroup.mat
% load the pig individual (originte from) for each cell
load cell_pigIndividual.mat

% remove cells with idx < 1
removeIdx = find(idx < 1);
clusterIdentifiers(removeIdx) = [];
idx(removeIdx) = [];
idxTxt(removeIdx) = [];
expressMat(:, removeIdx) = [];
coordinate(removeIdx, :) = [];
sampleID(removeIdx) = [];
originalSampleID(removeIdx) = [];
%% important: if reruning step 1, examine the localization of key marker genes:

% ensure that almost cells express endothelial cell markers (Manuscript Fig 1B-C)
threshold = 1;
plotMarkerOnUmap('PECAM1', threshold, gene, expressMat, coordinate);
plotMarkerOnUmap('KDR', threshold, gene, expressMat, coordinate);

% localization of vascular enthothelial cell markers (Manuscript Fig 1E-F)
threshold = 2;
plotMarkerOnUmap('CD34', threshold, gene, expressMat, coordinate);
plotMarkerOnUmap('PLVAP', threshold, gene, expressMat, coordinate);

% localization of lymphatic enthothelial cell markers (Manuscript Fig 1G-H)
threshold = 5;
plotMarkerOnUmap('CCL21', threshold, gene, expressMat, coordinate);
plotMarkerOnUmap('PROX1', threshold, gene, expressMat, coordinate);

% localization of 5 cell-cycle markers (Manuscript Fig 2B-F)
threshold = 0.5;
plotMarkerOnUmap('AURKB', threshold, gene, expressMat, coordinate);
threshold = 1;
plotMarkerOnUmap('ENSSSCG00000026302', threshold, gene, expressMat, coordinate); title('MKI67');

%% assign cluster number and name of the cluster. If reruning step 1, the cluster number assignment may change after examining the localization of key markers above
idx = 2*ones(size(clusterIdentifiers));
idx(clusterIdentifiers<1) = -1;
idxTxt = cell(size(idx));
idxTxt(1:length(idxTxt)) = {'NA'}; % cells with cluster -1, also called NA, don't belong to any cluster

idx( ismember( clusterIdentifiers, [19] ) ) = 1; 
idx( ismember( clusterIdentifiers, [25 26] ) ) = 3;
idx( ismember( clusterIdentifiers, [24] ) ) = 4;
idx( coordinate(:, 1) > 5 & clusterIdentifiers ~= 24) = 5;

idxTxt(idx == 1) = {'VEC1'};
idxTxt(idx == 2) = {'VEC2'};
idxTxt(idx == 3) = {'VEC3'};
idxTxt(idx == 4) = {'LEC1'};
idxTxt(idx == 5) = {'LEC2'};

% visualize the endothelial cell cluster (Manuscript Fig 1A)
figure, gscatter(coordinate(:, 1), coordinate(:, 2), idx); % manually delete '-1'
xlabel('umap 1'); ylabel('umap 2');
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off', 'FontSize', 16)

figure, gscatter(coordinate(:, 1), coordinate(:, 2), idxTxt); % manually delete 'NA'
xlabel('umap 1'); ylabel('umap 2');
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off', 'FontSize', 16)

% visualize cells by the pig groups and animals - negligible batch effect
% confirm
% by pig group (Supplemental Figure 2H)
figure, gscatter(coordinate(:, 1), coordinate(:, 2), sampleID ); h = legend('Location','eastoutside');
xlabel('umap 1'); ylabel('umap 2');
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off', 'FontSize', 16)
% by pig individual
figure, gscatter(coordinate(:, 1), coordinate(:, 2), originalSampleID ); h = legend('Location','eastoutside');
xlabel('umap 1'); ylabel('umap 2');
set(gca, 'linewidth', 1.5, 'XColor', 'k', 'Ycolor', 'k', 'TickDir', 'out', 'Box', 'off', 'FontSize', 16)

%% compute the statistics for each gene in each cluster
meanExp = zeros(length(gene), max(idx));
percenExp = zeros(length(gene), max(idx));
pFisher = ones(length(gene), max(idx));
oddRatio = ones(length(gene), max(idx));
meanExpGeneral = zeros(length(gene), 1);
for i = 1 : length(gene)
    meanExpGeneral(i) = mean( expressMat(i, :) );
    for j = 1 : max(idx)
        clusterIndex = find(idx == j);
        meanExp(i, j) = mean( expressMat(i, clusterIndex) );
        percenExp(i, j) = length( find( expressMat(i, clusterIndex) > 0 ) ) / length (clusterIndex);
        
        count1 = length ( find( expressMat(i, setdiff(1:length(idx), clusterIndex)) == 0 ) );
        count2 = length ( find( expressMat(i, clusterIndex) == 0 ) );
        count3 = length ( find( expressMat(i, setdiff(1:length(idx), clusterIndex)) > 0 ) );
        count4 = length ( find( expressMat(i, clusterIndex) > 0 ) );
        countNum = [ count1, count2; count3, count4];
        
        [h,pFisher(i, j),stats] = fishertest(countNum,'Tail','right','Alpha',0.01);
        oddRatio(i, j) = stats.OddsRatio;
    end
end

foldChange = zeros(length(gene), max(idx));
for clusterID = 1:max(idx)
    clusterIdx = find(idx == clusterID);
    nonClusterIdx = find(idx ~= clusterID);

    foldChange(:, clusterID) = mean( expressMat(:, clusterIdx)' )' ./ mean( expressMat(:, nonClusterIdx)' )';
end

% save the statistics in file Stat_EC_GOCycle_Only_NormalizedExpression
save ('Stat_EC_GOCycle_Only_NormalizedExpression.mat', 'idx', 'pFisher', 'meanExp', 'percenExp', 'oddRatio', 'meanExpGeneral', 'foldChange');

% save all results into a workspace file for Step 3
%save('EC_workspace.mat');


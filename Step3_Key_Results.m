%% Manuscript: Nguyen et al. Single-cell RNA sequencing analysis identifies 
% one subpopulation of endothelial cells that proliferates and another that 
% undergoes the endothelial-mesenchymal transition in regenerating pig
% hearts. Frontiers in Bioengineering and Biotechnology Tissue Engineering 
% and Regenerative Medicine. 2023. DOI: 10.3389/fbioe.2023.1257669. This
% manuscript examine endothelial cells and angiogenesis in pig hearts that
% recover from myocardial infarction induced 4 weeks after birth.

% This is Step 3 of the endothelial cells single nuclei RNA sequencing (snRNAseq)
% analyusis presented in the manuscript above. This step reproduce the key analytic figures.

% Results from 'Step 1' and 'Step 2' are already precomputed, so Step 3 can be executed
% without reruning Step 1. If rerunning step 1, be careful in cluster
% number assignment (see below).

%% clear the workspace (if needed)
clear all
close all
clc

%% load the endothelial cell workspace computed in Step 2
load EC_workspace.mat

%% the proportion of each endothelial cell clusters in each animal group
% all cluster proportions (Figure 1J)
PlotAllClusterProportion;

% cluster VEC1 & VEC3 (Figure 2A, Figure 4A)
PlotVEC_1_3_Proportion;

%% cluster heatmaps
% cell cycle & proliferation marker heatmap (Manuscript Figure 2G)
PlotCellCycleHeatmap

% endo-mesenchymal transition (EndMT) marker heatmap (simpler version of Manuscript Figure 4B)
PlotEndMTHeatmap

%% some helpful plot

% plot the expression / localization of cell expressing specific genes in
% 2D umap. The example below plots AURKB. To plot other genes, change the
% gene name in line 44. Change the expression threshold ( cells with
% expression <= threshold are represented in gray dots) in line 45
marker = 'AURKB';
threshold = 0.5;
plotMarkerOnUmap(marker, threshold, gene, expressMat, coordinate);

% bar plot of gene expression, group by cluster ID
plotMarkerExpressionBar(marker, gene, expressMat, idxTxt); title(marker); % RNA expression
% violin plot of gene expression, group by cluster ID
plotMarkerExpressionViolin(marker, gene, expressMat, idxTxt); title(marker); % RNA expression

% bar plot of gene expression, by animal group
plotMarkerExpressionBar(marker, gene, expressMat, sampleID); title(marker); % RNA expression
% violin plot of gene expression, by animal group
plotMarkerExpressionViolin(marker, gene, expressMat, sampleID); title(marker); % RNA expression

%% marker table for each cluster (manuscript supplemental table)
% change the cluster number in line 66, see follow 
% clusterID = 1: cluster VEC1
% clusterID = 2: cluster VEC2
% clusterID = 3: cluster VEC3
% clusterID = 4: cluster LEC1
% clusterID = 5: cluster LEC2

clusterID = 1;
% select cluster marker genes by the percentage of cells expressing the
% gene (percenExp), fold-change (foldChange) and p-value (pFisher)
% thresholds
markerIndex = find( percenExp(:, clusterID) > 0.3 & foldChange(:, clusterID) > 2 & pFisher(:, clusterID) < 0.05 );

% see the statisitcs in markerTable
markerTable =  table(gene(markerIndex), percenExp(markerIndex, clusterID), foldChange(markerIndex, clusterID), pFisher(markerIndex, clusterID)); 
markerTable.Properties.VariableNames = {'Gene', '%Cell expressing', 'Fold-change', 'p-value'};


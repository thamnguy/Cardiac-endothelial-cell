%% Manuscript: Nguyen et al. Single-cell RNA sequencing analysis identifies 
% one subpopulation of endothelial cells that proliferates and another that 
% undergoes the endothelial-mesenchymal transition in regenerating pig
% hearts. Frontiers in Bioengineering and Biotechnology Tissue Engineering 
% and Regenerative Medicine. 2023. DOI: 10.3389/fbioe.2023.1257669. This
% manuscript examine endothelial cells and angiogenesis in pig hearts that
% recover from myocardial infarction induced 4 weeks after birth.

% This file is Step 1 of analyzing the single-nuclei RNA sequencing (snRNAseq)
% data presented in the manuscript above. After this step, an autoencoder 
% and UMAP 2D visualization, which are needed for endothelial cells
% clustering, are computed.


% 'Step 0' of the analysis is the isolation of endothelial cell from the
% whole-heart snRNAseq data. This 'step 0' is presented in Nguyen et. al. 
% Cardiomyocyte cell-cycle regulation in neonatal large mammals: Single 
% Nucleus RNA-sequencing Data analysis via an Artificial-intelligence–based 
% pipeline. Frontiers in Bioengineering and Biotechnology Tissue Engineering 
% and Regenerative Medicine. 2022. https://doi.org/10.3389/fbioe.2022.914450

% Important: all files & folders, including subfiles and subfolders, inside
% the 'umapFileExchange' folder, must be added into Matlab via 'Set Path'.
% See
% https://www.mathworks.com/matlabcentral/answers/247180-how-can-i-add-all-subfolders-to-my-matlab-path#accepted_answer_194998.
% for instruction
% The folder is an earlier version of Matlab UMAP toolkit: 
% https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap#version_history_tab

% Result of this steps is already stored; therefore, to visualize the
% results showed in the manuscript above, this step is not needed. 

%% clear the workspace
clear all
close all
clc

%% a. load the endothelial cell expression data
% load the endothelial cell snRNAseq data. The data is already normalized in 'step 0' 
load EC_expressMat.mat 

% load the entire gene list in the data
load allGene.mat;

% load the cell-cycle-specific gene list
load cellCycleGene.mat;

%% compute the Autoencoder model that would cluster the endothelial cells
cellCycleIndex = find( ismember(gene, cellCycleGene)==1 );
autoencodeExpressMat = full( log2( expressMat( cellCycleIndex , : ) + 1) );

autoenc = trainAutoencoder(autoencodeExpressMat, 'UseGPU', true); 
% if without GPU or low-memory GPU, run
% autoenc = trainAutoencoder(autoencodeExpressMat, 'UseGPU', false); 

%% cluster the endothelial cells
Z = encode(autoenc,expressMat);

% Matlab UMAP toolkit: 
% https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap#version_history_tab
[coordinate, umap, clusterIdentifiers]=run_umap(Z', 'n_epoch', 500);
clusterIdentifiers = clusterIdentifiers';

% important: different rerun of the Autoencoder and UMAP may lead to
% different cluster number assignment. Therefore, Step 1 stops at this
% point. Final cluster assignment should be done in Step 2 with careful
% examinination of cluster number and marker genes expression.

%% save the result
save allEC_GOCycle_Autoencoder.mat autoenc -mat; % the autoencoder model
save allEC_GOCycle_Endocde.mat Z -mat; % the autoencoder processing (embedding) of endothelial cells
save allEC_GOCycle_Umap.mat coordinate -mat; % all cells 2D UMAP for visualization
save allEC_GOCycle_UmapClusterID.mat clusterIdentifiers -mat; %pre-clustering
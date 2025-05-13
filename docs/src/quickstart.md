```@meta
CurrentModule = SGCRNAs
```

```@index
```

# Quick Start
SGCRNAs consists of five functions:
- CGM: Performs preprocessing for clustering.
- SpectralClustering: Carries out the clustering.
- SetNetwork: Prepares data for network graph visualisation.
- DrawNetwork: Generates the network graph.

## Example
```
using CSV, DataFrames
using SGCRNAs

# Load Data
Data = CSV.read("Result/Norm/normalizedCounts_coding.tsv", header=1, comment="#", delim='\t', DataFrame);

# Pre-proccessing
CorData, GradData = CGM(Data.Symbol, Matrix(Data[:,5:end]), fn="Result/coding-FTEST", mode=:FTEST);

# Clustering
clust, pos, edge_data = SpectralClustering(CorData, GradData);

# Draw Network
d = 1; k = maximum(clust[d]);
nw, new_pos, cnctdf, new_clust, score = SetNetwork(edge_data, clust[d], pos, il=collect(1:k));
DrawNetwork("Result/Fig/AllNetWork-0.5.png", nw, new_pos, cnctdf, new_clust, k, node_scores=score, edge_mode=:ALL, edge_threshold=0.5);

# module-phenomenon correlation
Phen = CSV.read("SraRunTable.csv", header=1, comment="#", delim=',', DataFrame);
sort!(Phen);
Data = innerjoin(Data, DataFrame(Symbol=names(edge_data)), on=:Symbol, order=:right);
SGCRNAs.CorPhenMod(Data[:,5:end], Phen[:,[2,3,5]], new_clust, "Result/Fig/CorPhenMod.png");
```
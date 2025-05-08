# SGCRNAs [![Build Status](https://github.com/C37H41N2O6/SGCRNAs/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/C37H41N2O6/SGCRNAs/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage Status](https://coveralls.io/repos/github/C37H41N2O6/SGCRNAs/badge.svg?branch=main)](https://coveralls.io/github/C37H41N2O6/SGCRNAs?branch=main) [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://C37H41N2O6.github.io/SGCRNAs.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://C37H41N2O6.github.io/SGCRNAs.jl/dev/)

## Installation
Install the package using Julia's package manager:
```julia
using Pkg; Pkg.add("SGCRNAs")
```

## Usage
Here's a basic example:
```julia
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

## Documentation
Full documentation is available at: https://C37H41N2O6.github.io/SGCRNAs.jl

## License
This project is licensed under the MIT License. See LICENSE for details.

## Related Publication
If you use this package in your research, please cite the following paper:
Author(s), Title of the Paper, Journal/Conference, Year.
[DOI or URL to the publication]
```
@article{YourCitationKey,
  author    = {Author, First and Author, Second},
  title     = {Title of the Paper},
  journal   = {Journal Name},
  year      = {202X},
  volume    = {X},
  number    = {Y},
  pages     = {Z-ZZ},
  doi       = {10.xxxx/xxxxx}
}
```

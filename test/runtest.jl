# test/runtest.jl

using CSV, DataFrames
using Random
using SGCRNAs
using Test

num_genes=100; num_samples=20; num_modules=4;
Random.seed!(42)
genes_per_module = div(num_genes, num_modules)
base_profiles = [randn(num_samples) for _ in 1:num_modules]
expression_data = Matrix{Float64}(undef, num_genes, num_samples)
for i in 1:num_genes
    clust = div(i - 1, genes_per_module) + 1
    noise = 0.2 * randn(num_samples)
    expression_data[i, :] = base_profiles[clust] .+ noise
end
gene_names = ["Gene_$(i)" for i in 1:num_genes]

@testset "Function tests" begin
    @test_nowarn begin
        CorData, GradData = CGM(gene_names, expression_data);
    end
    @test_nowarn begin
        clust, pos, edge_data = SpectralClustering(CorData, GradData);
    end
    d = 1; k = maximum(clust[d]);
    @test_nowarn begin
        nw, new_pos, cnctdf, new_clust, score = SetNetwork(edge_data, clust[d], pos, il=collect(1:k));
    end
    @test_nowarn begin
        DrawNetwork("./AllNetWork-0.5.png", nw, new_pos, cnctdf, new_clust, k, node_scores=score, edge_mode=:ALL, edge_threshold=0.5)
    end
    @test_nowarn begin
        CorPhenMod(DataFrame(expression_data, :auto), Matrix(expression_data')[:,[3,1,4]], new_clust, "./CorPhenMod.png")
    end
end

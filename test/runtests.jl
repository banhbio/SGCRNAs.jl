using CSV, DataFrames
using LinearAlgebra, Distributions
using Random
using SGCRNAs
using Test

num_genes=500; num_samples=20; num_modules=4;
Random.seed!(42)
genes_per_module = div(num_genes, num_modules)
# common_trend = randn(num_samples)
# base_profiles = [common_trend .* rand() .+ 0.1 * randn(num_samples) for _ in 1:num_modules]
# expression_data = Matrix{Float64}(undef, num_genes, num_samples)
# for i in 1:num_genes
#     clust = div(i - 1, genes_per_module) + 1
#     noise = 0.01 * randn(num_samples)
#     expression_data[i, :] = base_profiles[clust] .+ noise
# end
gene_names = ["Gene_$(i)" for i in 1:num_genes]
ρ = 0.8
Σ = [ρ^abs(i-j) for i in 1:num_samples, j in 1:num_samples]
dist = MvNormal(zeros(num_samples), Symmetric(Σ))
expression_data = [rand(dist) for _ in 1:num_genes]
expression_data = reduce(vcat, expression_data)

@testset "Function tests" begin
    CorData, GradData = (0, 0);
    @test_nowarn begin
        CorData, GradData = CGM(gene_names, expression_data);
    end
    clust, pos, edge_data = (0, 0, 0);
    @test_nowarn begin
        clust, pos, edge_data = SpectralClustering(CorData, GradData);
    end
    d = 1; k = maximum(clust[d]);
    nw, new_pos, cnctdf, new_clust, score = (0, 0, 0, 0, 0);
    @test_nowarn begin
        nw, new_pos, cnctdf, new_clust, score = SetNetwork(edge_data, clust[d], pos, il=collect(1:k));
    end
    @test_nowarn begin
        DrawNetwork("./AllNetWork-0.5.png", nw, new_pos, cnctdf, new_clust, k, node_scores=score, edge_mode=:ALL, edge_threshold=0.5)
    end
    @test_nowarn begin
        CorPhenMod(DataFrame(expression_data, :auto), DataFrame(Matrix(expression_data')[:,[3,1,4]], :auto), new_clust, "./CorPhenMod.png")
    end
end

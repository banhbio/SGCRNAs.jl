module SGCRNAs
    using Dates
    using CSV, DataFrames
    using KrylovKit
    using StatsBase
    using LinearAlgebra, Statistics, MultivariateStats, Distributions, KernelDensity
    using ParallelKMeans, Clustering
    using Random, UMAP
    using Graphs, Colors


    ##### correlation & gradient matrix calculation #####
        """
        # arguments
        - gene::Vector: gene name list
        - data::Matrix: gene expression matrix
        - threshold::Float64: value for remove genes with more than a certain number of zeros; default: 0.5
        - mode::Symbol: mode of measurement errors elimination
          - :NONE -> measurement error is not considered (Select this option when there are sufficient number of samples; Default)
          - :LESS -> defined as the value below the mode
          - :SIGMA -> defined as the value below the 2σ(mode is considered as σ)
          - :FTEST -> defined as the value below the significantly different from the mode by pval
        - binSize::Float64: histogram bin size used to determine measurement error when mode is other than :NONE; default: 0.01
        - pval::Float64: p-value for determining measurement error using :FTEST & Statistical tests of correlation coefficients; default: 0.05
        - power::Float64: Power in statistical tests of correlation coefficients; default: 0.8
        # returns
        - CorData: correlation matrix
        - GradData: gradient matrix
        """
        function CGM(gene::Vector, data::Matrix; threshold::Float64=0.5, mode::Symbol=:NONE, binSize::Float64=0.01, pval::Float64=0.05, power::Float64=0.8)
            SmplNum = size(data, 2)
        
            # Remove genes with more than a certain number of zeros
            Q = map(x -> sum(x .== 0.0), eachrow(data))
            Gene = gene[Q .< SmplNum * threshold]
            Data = data[Q .< SmplNum * threshold, :]
        
            # Eliminates measurement errors
            if mode != :NONE
                # Get the mode of the coefficient of variation
                cvList = map(x -> std(x, corrected=false), eachrow(Data)) ./ map(mean, eachrow(Data))
                freqCurve = pdf(kde(cvList), [0.0:binSize:ceil(maximum(cvList))])[1]
                cvMode = binSize * (argmax(freqCurve) - 1)

                if mode == :LESS
                    # Remove cvMode and below
                elseif mode == :SIGMA
                    # Assume cvMode as σ and remove 2σ
                    cvMode *= 2
                elseif mode == :FTEST
                    # Delete cv that is below the significance level in F test
                    cvMode *= sqrt(quantile(FDist(SmplNum-1, SmplNum-1), 1 - pval))
                end
        
                Gene = Gene[cvList .> cvMode]
                Data = Data[cvList .> cvMode, :]
            end
        
            # Calculate covariance and standard deviation
            Avg = map(mean, eachrow(Data))
            Dist = Data .- Avg
            Covar = (Dist * Dist')
            Std = map(x -> sqrt(sum(x .^ 2)), eachrow(Dist))
            # Calculate correlation coefficients in batches
            CorMat = Covar ./ (Std * Std')
            # Calculate gradients in batches
            VarVec = map(x -> sum(x .^ 2), eachrow(Dist))
            GeneNum = length(VarVec)
            VarMat = zeros(GeneNum, GeneNum)
            for i in 1:GeneNum
                VarMat[i, 1:i] .= VarVec[i]
            end
            VarMat += VarMat'
            for i in 1:GeneNum
                VarMat[i, i] -= VarVec[i]
            end
            CorTerm = Avg' ./ Avg
            CorTerm = triu(CorTerm) + triu(CorTerm, 1)'
            # CorTerm = (CorTerm .+ CorTerm') ./ 2
            GradMat = (Covar ./ VarMat) .* CorTerm

            # Set all but statistically significant correlation coefficients to zero
            d = Normal()
            Za = quantile(Normal(), 1-pval/2)
            Zb = quantile(Normal(), power)
            z = exp(2*(Za+Zb) / sqrt(SmplNum-3))
            r = (z-1) / (z+1)
            CorMat .*= (abs.(CorMat) .>= r)

            # Conversion to data frame
            CorData = hcat(
                        DataFrame(Symbol=Gene),
                        DataFrame(CorMat, Gene)
                    )
            CorData = coalesce.(CorData, 0.0)
            replace!.(eachcol(CorData), NaN => 0.0)
            GradData = hcat(
                        DataFrame(Symbol=Gene),
                        DataFrame(GradMat, Gene)
                    )
            GradData = coalesce.(GradData, 0.0)
            replace!.(eachcol(GradData), NaN => 0.0)
                
            # Genes that did not correlate with any of the genes were removed
            CorData = CorData[sum.(eachrow(CorData[:,2:end])) .!= 0.0, :]
            GradData = innerjoin(CorData[:,[:Symbol]], GradData, on=:Symbol)

            # alignment
            sort!(CorData, :Symbol)
            CorData = CorData[:, sort(CorData.Symbol)]
            sort!(GradData, :Symbol)
            GradData = GradData[:, sort(GradData.Symbol)]
        
            # Keep 0 between selves.
            for i in 1:ncol(CorData)
                CorData[i, i] = 0.0
                GradData[i, i] = 0.0
            end

            return CorData, GradData
        end
        export CGM
    ##### correlation matrix calculation #####

    ##### Laplacian matrix calculation #####
        function LaplacianMatrix(mat::Matrix, normFlg::Bool, randNormFlg::Bool)
            nodeScores = sum.(eachrow(mat))
            matD = diagm(nodeScores)
            matL = matD .- mat
            if(normFlg)
                matL = 1.0I(size(mat,2)) .- sqrt(inv(matD)) * matL * sqrt(inv(matD))
            elseif(randNormFlg)
                matL = 1.0I(size(mat,2)) .- inv(matD) * matL
            end

            return matL
        end
    ##### Laplacian matrix calculation #####

    ##### clustering #####
        function Clustering_Eigen(matL::Matrix, maxK::Int64, normFlg::Bool)
            eigVals, eigVecs, eigInfo = eigsolve(matL, maxK+10, :SR, krylovdim=5*maxK)
            eigVals = Real.(eigVals)
            eigVecs = Real.(reduce(hcat, eigVecs)')

            # calculate normalized gap
            normGaps = [(eigVals[k+1] - eigVals[k]) / eigVals[k] for k in 1:(length(eigVals)-1)]
            sortedGaps = sortperm(normGaps, rev=true)
            k = sortedGaps[1] == 1 ? sortedGaps[2] : sortedGaps[1]

            embedding = eigVecs[1:k, :]
            if (normFlg)
                buf = map(x -> x ./ sum(x .^ 2), eachrow(embedding))
                embedding = reduce(hcat, buf)'
            end
            
            return Matrix(embedding)
        end
    ##### clustering #####

    ##### SpectralClustering #####
        function Clustering_Main(df::DataFrame, itr::Int64, seed::Int64, pcas::Int64, normFlg::Bool, randNormFlg::Bool)
            RndSeed = Random.seed!(seed)

            matL = LaplacianMatrix(Matrix(df), normFlg, randNormFlg)
            emb = Clustering_Eigen(matL, pcas, normFlg)
            k = size(emb,1)

            # k-means clustering by automatic k-value determination
            if size(emb, 2) < 100
                res = ParallelKMeans.kmeans(Hamerly(), emb, k, max_iters=itr, rng=RndSeed)
            elseif size(emb, 2) > 10000
                res = ParallelKMeans.kmeans(Elkan(), emb, k, max_iters=itr, rng=RndSeed)
            else
                res = ParallelKMeans.kmeans(Yinyang(), emb, k, max_iters=itr, rng=RndSeed)
            end

            return emb, res.assignments
        end
        """
        # arguments
        - cor::DataFrame: dataframe of correlation matrix (return value of CGM())
        - grad::DataFrame: dataframe of gradient matrix (return value of CGM())
        - tNodeNum::Int64: threshold of sub-cluster node number; default: 100
        - depthMaxv: Depth of sub-clusters; default: 5
        - pcas::Int64: pca dimention; default: 99
        - itr::Int64: number of trials; default: 300
        - seed::Int64: seed value of random number; default: 42 (Answer to the Ultimate Question of Life, the Universe, and Everything)
        - nNeighbors::Int64: UMAP parameter; default: 40
        - minDist::Float64: UMAP parameter; default: 0.1
        - normFlg::Bool: Whether to symmetrically normalize the Laplacian matrix; default: true
        - randNormFlg::Bool: Whether to random walk normalize the Laplacian matrix; default: false
        # returns
        - clust: cluster number of each gene
        - pos: gene position for drawing network
        - edgeScore: edge score for drawing network
        """
        function SpectralClustering(cor::DataFrame, grad::DataFrame; tNodeNum::Int64=100, depthMax::Int64=5, pcas::Int64=99, itr::Int64=300, seed::Int64=42, nNeighbors::Int64=40, minDist::Float64=0.1, normFlg::Bool=true, randNormFlg::Bool=false)
            rowNum = size(cor, 1)
            df = ((1 .+ cor) ./ 2) .* exp.(-1 .* abs.(log.(abs.(grad))))
            # Laplacian matrix calculation
            emb, clust = Clustering_Main(df, itr, seed, pcas, normFlg, randNormFlg)
            clustData = [clust]
            kMax = maximum(clust); d = 0;
            while ((maximum(map(x -> sum(clustData[d+1] .== x), 1:kMax)) > tNodeNum) & (d < depthMax))
                append!(clustData, deepcopy([clustData[d+1]]))
                for k in 1:kMax
                    Q = (clustData[d+2] .== k)
                    subDF = df[Q, Q]
                    if size(subDF, 2) > tNodeNum
                        _, subClust = Clustering_Main(subDF, itr, seed, pcas, normFlg, randNormFlg)
                        subClust .-= 1
                        subClust[subClust .!= 0] .+= kMax
                        subClust[subClust .== 0] .= k
                        clustData[d+2][Q] = subClust
                        kMax = maximum(subClust)
                    end
                end
                d += 1
            end

            if (size(emb, 1) > 2)
                embedding = umap(emb, 2; n_neighbors=nNeighbors, min_dist=minDist)
            else
                embedding = emb
            end

            return clustData, Matrix(permutedims(embedding)), cor .* exp.(-1 .* abs.(log.(abs.(grad))))
        end
        export SpectralClustering
    ##### SpectralClustering #####

    ##### draw network #####
        """
        # arguments
        - df::DataFrame: dataframe of correlation matrix (return value of CGM())
        - clust::Vector{Int64}: cluster number of each gene (one of return value of SpectralClustering())
        - pos::Matrix: gene position for drawing network (one of return value of SpectralClustering())
        - il::Vector: module number list which you want to draw
        # returns
        - nw: undirected graph
        - pos: node position
        - cnctdf: converted correlation matrix
        - clust: cluster number of each gene in network
        - score: node scores
        """    
        function SetNetwork(df::DataFrame, clust::Vector{Int64}, pos::Matrix; il::Vector=[])
            ##### preliminaries #####
                # All clusters you want to draw if none are specified.
                if length(il) == 0
                    il = sort(unique(clust))
                end
                # Extract genes present in the cluster you want to draw
                Q1 = (clust .== il[1])
                if length(il) > 1
                    for i in 2:length(il)
                        Q1 .|= (clust .== il[i])
                    end
                end
                cnctdf = deepcopy(Matrix(df[Q1, Q1]))
                gene_list = names(df)[Q1]
                Q2 = (map(sum, eachrow(cnctdf)) .!= 0.0)
                cnctdf = cnctdf[Q2, Q2]
                gene_list = gene_list[Q2]
                gene_num = length(gene_list)

                # Convert to upper triangular matrix
                triu!(cnctdf)
                cnctdf = DataFrame(hcat(gene_list,cnctdf), vcat(["Symbol"],gene_list))
                cnctdf = stack(cnctdf, 2:ncol(cnctdf))
                rename!(cnctdf, [:e1,:e2,:cor])
                # Remove correlations in the same gene and duplicate combinations
                ## Correlation coefficients between themselves are set to 0.
                ## The overlapping combinations have a correlation coefficient of zero
                ## due to the conversion to an upper triangular matrix.
                cnctdf = cnctdf[cnctdf.cor .!= 0.0, :]
                sort!(cnctdf, :e1)
            
                # Assign node numbers to genes
                buf = DataFrame(e1=sort(unique(vcat(cnctdf.e1, cnctdf.e2))))
                buf[!, :i1] = collect(1:nrow(buf))
                cnctdf = innerjoin(cnctdf, buf, on=:e1)
                rename!(buf, [:e2,:i2])
                cnctdf = innerjoin(cnctdf, buf, on=:e2)

                # Assign module numbers to genes
                buf = DataFrame(e1=gene_list, m1=clust[Q1][Q2])
                cnctdf = innerjoin(cnctdf, buf, on=:e1)
                rename!(buf, [:e2,:m2])
                cnctdf = innerjoin(cnctdf, buf, on=:e2)
            ##### preliminaries #####

            ##### Graph Generation #####
                # undirected graph
                nw = SimpleGraph(gene_num)
                for i in 1:nrow(cnctdf)
                    add_edge!(nw, cnctdf.i1[i], cnctdf.i2[i])
                end
                # weighted degree by edge value
                score = vec(sum(abs.(Matrix(df[Q1, Q1][Q2, Q2])), dims=2))
            ##### Graph Generation #####

            return nw, pos[Q1, :][Q2, :], cnctdf, clust[Q1][Q2], score
        end
        export SetNetwork
        """
        # arguments
        - fn: figure save name
        - nw: network graph (one of return value of SetNetwork())
        - pos: node position (one of return value of SetNetwork())
        - cnctdf: converted correlation matrix (one of return value of SetNetwork())
        - clust: cluster number of each gene in network (one of return value of SetNetwork())
        - k: number of clusters
        - node_scores: weight of node
        - node_labels: label of node
        - node_scaler: multiple for node diameter adjustment; Default: 100
        - edge_mode: mode of edges to be drawn
          - :ALL -> All edges are drawn (Default)
          - :N -> Only draw edges with negative values
          - :P -> Only draw edges with positive values
        - edge_threshold: Threshold value of edges to be drawn; Default: 0.0
        - edge_scaler: multiple for edge thickness adjustment; Default: 5
        - x_size, y_size: Size of the drawing area; Default: 100, 100
        """
        function DrawNetwork(fn::String, nw::SimpleGraph, pos::Matrix, cnctdf::DataFrame, clust::Vector{Int64}, k::Int64; node_scores::Vector{}=[], node_labels::Vector{}=[], node_color::Vector=[], node_scaler::Int64=100, edge_mode::Symbol=:ALL, edge_threshold::Float64=0.0, edge_scaler::Int64=5, x_size::Int64=50, y_size::Int64=50)
            gene_num = nv(nw)
        
            # node設計
            if length(node_scores) == 0
                node_scores = repeat([1], gene_num)
            end
            node_sizes = node_scores .* node_scaler
            if length(node_color) == 0
                clust_num = length(unique(vcat(cnctdf.m1, cnctdf.m2)))
                color_list = range(LCHuv(65,100,15), stop=LCHuv(65,100,375), length=k+1)

                node_color = [color_list[i] for i in clust]
            end
            if length(node_labels) == 0
                node_labels = repeat([""], gene_num)
            end
            
            # edge設計
            edge_colors = RGBA.(1.0, 0.3, 0.0, cnctdf.cor)
            edge_colors[cnctdf.cor .< 0.0] = RGBA.(0.0, 0.35, 1.0, abs.(cnctdf.cor[cnctdf.cor .< 0.0]))
            Q = (abs.(cnctdf.cor) .< edge_threshold)
            if edge_mode == :P
                Q = (cnctdf.cor .< 1*edge_threshold)
            elseif edge_mode == :N
                Q = (cnctdf.cor .> -1*edge_threshold)
            end
            if sum(Q) == nrow(cnctdf)
                println("NoEdge")
                return nothing
            else
                edge_colors[Q] .= RGBA(0.0, 0.0, 0.0, 0.0)
            end
            edge_sizes = edge_scaler .* abs.(cnctdf.cor)

            fig = gplot(
                        nw, pos[:, 1] .+ minimum(pos[:, 1]), pos[:, 2] .+ minimum(pos[:, 2]),
                        nodesize=node_sizes, nodefillc=node_color, nodelabel=node_labels,
                        edgestrokec=edge_colors, edgelinewidth=edge_sizes
                    )
            Compose.draw(PNG(fn, x_size*cm, y_size*cm), fig)

            return nothing
        end
        export DrawNetwork
    ##### draw network #####

    ##### Correlation of Phenomenon and Modules #####
        """
        # arguments
        - df1::DataFrame: dataframe of gene expression
        - df2::DataFrame: dataframe of Phenomenon
        - clust::Vector{Int64}: cluster number of each gene (one of return value of SpectralClustering())
        - fn::String: fig save name
        - cor_mode::Symbol: mode of caluclation of correlation coefficient
          - :ALL -> All three types are drawn (default)
          - :A_AVG -> all gene average
          - :P_AVG -> positive correlation gene average
          - :N_AVG -> negative correlation gene average
        """
        function CorPhenMod(df1::DataFrame, df2::DataFrame, clust::Vector{Int64}, fn::String; cor_mode::Symbol=:ALL)
            kuni =  sort(unique(clust))
            knum = length(kuni)

            CorList = []
            for i in 1:ncol(df2)
                push!(CorList, [[] for i=1:knum])
            end
            for k in 1:knum
                for i in 1:ncol(df2)
                    buf = df1[clust .== kuni[k], :]
		            Result_each = zeros(nrow(buf))
                    for j in 1:nrow(buf)
                        Result_each[j] = cor(Array(buf[j,:]), df2[:,i])
                    end
		            CorList[i][k] = Result_each
                end
            end

            x = collect(1:ncol(df2))
            y = collect(knum:-1:1)

        end
        export CorPhenMod
    ##### Correlation of Phenomenon and Modules #####
end

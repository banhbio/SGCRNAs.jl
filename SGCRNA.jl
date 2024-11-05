module SGCRNAs
    using Dates
    using CSV, DataFrames
    using KrylovKit
    using StatsBase
    using LinearAlgebra, Statistics, MultivariateStats, Distributions, KernelDensity
    using ParallelKMeans, Clustering
    using Random, UMAP
    using RCall
    using NetworkLayout, GraphPlot, Graphs, Colors
    using Compose, CairoMakie, Fontconfig
    using JLD2

    ##### get rank of element #####
        function ElmRank(arr::Vector, revflg::Bool)
            sorted_indices = sortperm(arr, rev=revflg)
            sorted_arr = arr[sorted_indices]
            ranks = similar(arr)

            rank = 1
            for i in 1:length(arr)
                if i > 1 && sorted_arr[i] != sorted_arr[i - 1]
                    rank = i
                end
                ranks[sorted_indices[i]] = rank
            end

            return ranks
        end
    ##### get rank of element #####

    ##### correlation matrix calculation #####
        """
        argument
        gene: gene name list
        data: gene expression matrix
        fn: prefix of save file name. defalut: do not save
            output file name: fn_cor.tsv
        threshold: value for remove genes with more than a certain number of zeros
        mode: mode of measurement errors elimination
                :NONE -> measurement error is not considered (Select this option when there are sufficient number of samples; Default)
                :LESS -> defined as the value below the mode
                :SIGMA -> defined as the value below the 2σ(mode is considered as σ)
                :FTEST -> defined as the value below the significantly different from the mode by pval
        binSize: histogram bin size used to determine measurement error when mode is other than :NONE
        pval: p-value for determining measurement error using :FTEST & Statistical tests of correlation coefficients
        power: Power in statistical tests of correlation coefficients.

        returns
        CorData: correlation matrix
        GradData: gradient matrix
        """
        function CCM(gene::Vector, data::Matrix; fn::String="", threshold::Float64=0.5, mode::Symbol=:NONE, binSize::Float64=0.01, pval::Float64=0.05, power::Float64=0.8)
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

                # plot freqCurve
                xs = collect(0.0:binSize:ceil(maximum(cvList)))
                f = Figure(size = (1300, 400))
                ax1 = Axis(f[1, 1], title = ":LESS", ylabel="probability density", xautolimitmargin=(0.0f0, 0.0f0), yautolimitmargin=(0.0f0, 0.0f0))
                ax2 = Axis(f[1, 2], title = ":SIGMA", xautolimitmargin=(0.0f0, 0.0f0), yautolimitmargin=(0.0f0, 0.0f0), yticklabelsvisible=false)
                ax3 = Axis(f[1, 3], title = ":FTEST", xautolimitmargin=(0.0f0, 0.0f0), yautolimitmargin=(0.0f0, 0.0f0), yticklabelsvisible=false)
                linkyaxes!(ax1, ax2)
                linkyaxes!(ax1, ax3)
                hidespines!(ax1, :t, :r)
                hidespines!(ax2, :t, :r)
                hidespines!(ax3, :t, :r)
                hidedecorations!(ax1, label=false, ticklabels=false, ticks=false)
                hidedecorations!(ax2, label=false, ticklabels=false, ticks=false)
                hidedecorations!(ax3, label=false, ticklabels=false, ticks=false)
                barplot!(ax1, xs, freqCurve, color=[x <= cvMode ? "red" : "blue" for x in xs], gap=0)
                barplot!(ax2, xs, freqCurve, color=[x <= cvMode*2 ? "red" : "blue" for x in xs], gap=0)
                barplot!(ax3, xs, freqCurve, color=[x <= cvMode*sqrt(quantile(FDist(SmplNum-1, SmplNum-1), 1 - pval)) ? "red" : "blue" for x in xs], gap=0)
                Label(f[end+1, :], text="standard deviation")
                CairoMakie.save(fn * "_freqCurve.png", f)

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
            VarMat = Matrix{Float64}(undef, GeneNum, GeneNum)
            for i in 1:GeneNum
                VarMat[i, 1:i] .= VarVec[i]
            end
            VarMat += VarMat'
            for i in 1:GeneNum
                VarMat[i, i] -= VarVec[i]
            end
            GradMat = Covar ./ VarMat

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

            save_object(fn * "_cor.jld2", CorData)
            save_object(fn * "_grad.jld2", GradData)

            return CorData, GradData
        end
        export CCM
    ##### correlation matrix calculation #####

    ##### Laplacian matrix calculation #####
        function LaplacianMatrix(df::DataFrame, normFlg::Bool, randNormFlg::Bool)
            mat = Matrix(df)
            nodeScores = sum.(eachrow(mat))
            matD = diagm(nodeScores)
            matL = matD .- mat
            if(normFlg)
                matL = sqrt(inv(matD)) * matL * sqrt(inv(matD))
            elseif(randNormFlg)
                matL = inv(matD) * matL
            end

            return matL
        end
    ##### Laplacian matrix calculation #####

    ##### clustering #####
        function Clustering_Eigen(matL::Matrix, k::Int64, normFlg::Bool)
            eigVals, eigVecs, eigInfo = eigsolve(matL, k+10, :SR, krylovdim=5*k)
            eigVals = Real.(eigVals)
            eigVecs = Real.(reduce(hcat, eigVecs))
            eigVecs = eigVecs[:, 1e-8 .< eigVals]

            embedding = eigVecs'[1:k, :]
            if (normFlg)
                buf = map(x -> x ./ sum(x .^ 2), eachrow(embedding))
                embedding = reduce(hcat, buf)'
            end
            
            return Matrix(embedding)
        end
        function Clustering_Main(emb::Matrix, seed::Int64, k::Int64, itr::Int64)
            RndSeed = Random.seed!(seed)

            embedding = deepcopy(emb[1:k, :])
            # k-means clustering by automatic k-value determination
            if size(embedding, 2) < 100
                res = ParallelKMeans.kmeans(Hamerly(), embedding, k, max_iters=itr, rng=RndSeed)
            elseif size(embedding, 2) > 10000
                res = ParallelKMeans.kmeans(Elkan(), embedding, k, max_iters=itr, rng=RndSeed)
            else
                res = ParallelKMeans.kmeans(Yinyang(), embedding, k, max_iters=itr, rng=RndSeed)
            end
            

            score = [0.0,0.0,0.0,0.0]
            score[1] = clustering_quality(embedding, res.centers, res.assignments, quality_index=:xie_beni)
            score[2] = clustering_quality(embedding, res.centers, res.assignments, quality_index = :davies_bouldin)
            score[3] = clustering_quality(embedding, res.assignments, quality_index=:silhouettes)
            score[4] = clustering_quality(embedding, res.centers, res.assignments, quality_index = :calinski_harabasz)
            
            return embedding, res, score
        end
    ##### clustering #####

    ##### SpectralClustering #####
        function SCSub(emb::Matrix, kmin::Int64, kstep::Int64, kmax::Int64, itr::Int64, seed::Int64)
            resDf = DataFrame(k=[], assign=[], xbi=[], dbi=[], si=[], chi=[])
            for k in kmin:kstep:kmax
                embedding, res, score = Clustering_Main(emb, seed, k, itr)
                push!(resDf, hcat([k], [res.assignments], score'))
            end
            resDf.k = float.(resDf.k)
            resDf[!, 3:6] = float.(resDf[:, 3:6])

            # auto detect value
            RankList = [ElmRank(resDf.si,true),ElmRank(resDf.chi,true),ElmRank(resDf.xbi,false),ElmRank(resDf.dbi,false)]
            RankList = hcat(RankList...)
            Score = map(prod, eachrow(RankList))
            minK = argmin(Score)

            return resDf.assign[minK]
        end
        """
        argument
        cor: dataframe of correlation matrix (return value of CCM())
        grad: dataframe of gradient matrix (return value of CCM())
        tNodeNum: threshold of sub-cluster node number; default: 100
        depthMax: Depth of sub-clusters; default: 5
        kmin,kstep,kmax: test range of clusters number; default: 3, 1, 15
        itr: number of trials; default: 300
        seed: seed value of random number; default: 42 (Answer to the Ultimate Question of Life, the Universe, and Everything)
        nNeighbors: UMAP parameter; default: 40
        minDist: UMAP parameter; default: 0.01
        normFlg: Whether to symmetrically normalize the Laplacian matrix; default: false
        randNormFlg: Whether to random walk normalize the Laplacian matrix; default: false

        returns
        clust: cluster number of each gene
        pos: gene position for drawing network
        edgeScore: edge score for drawing network
        """
        function SpectralClustering(cor::DataFrame, grad::DataFrame; tNodeNum::Int64=100, depthMax::Int64=5, kmin::Int64=3, kstep::Int64=1, kmax::Int64=15, itr::Int64=300, seed::Int64=42, nNeighbors::Int64=40, minDist::Float64=0.01, normFlg::Bool=true, randNormFlg::Bool=false)
            rowNum = size(cor, 1)
            df = ((1 .+ cor) ./ 2) .* exp.(-1 .* abs.(log.(abs.(grad))))
            # Laplacian matrix calculation
            matL = LaplacianMatrix(df, normFlg, randNormFlg)
            emb = Clustering_Eigen(matL, kmax, normFlg)
            clust = SCSub(emb, kmin, kstep, kmax, itr, seed)
            clustData = [clust]
            kNum = maximum(clust[1]); d = 0;
            while ((maximum(map(x -> sum(clustData[d+1] .== x), 1:kNum)) > tNodeNum) & (d < depthMax))
                append!(clustData, deepcopy([clustData[d+1]]))
                for k in 1:kNum
                    subEmb = emb[:, clustData[d+2] .== k]
                    if size(subEmb, 2) > tNodeNum
                        subClust = SCSub(subEmb, kmin, kstep, kmax, itr, seed) .- 1
                        subClust[subClust .!= 0] .+= kNum
                        subClust[subClust .== 0] .= k
                        clustData[d+2][clustData[d+2] .== k] = subClust
                        kNum = maximum(subClust)
                    end
                end
                d += 1
            end

            kNum = maximum(clustData[1])
            embedding = umap(emb[1:kNum, :], 2; n_neighbors=nNeighbors, min_dist=minDist)

            # # plot corrected SSE
            # ks = collect(kmin:kstep:kmax)
            # f = Figure(size = (1125, 750))
            # ax1 = Axis(f[1, 1], title = "Xie-Beni")
            # ax2 = Axis(f[1, 2], title = "Davies-Bouldin")
            # ax3 = Axis(f[2, 1], title = "silhouettes")
            # ax4 = Axis(f[2, 2], title = "Calinski-Harabasz")
            # hidespines!(ax1, :t, :r)
            # hidespines!(ax2, :t, :r)
            # hidespines!(ax3, :t, :r)
            # hidespines!(ax4, :t, :r)
            # hidedecorations!(ax1, label=false, ticklabels=false, ticks=false)
            # hidedecorations!(ax2, label=false, ticklabels=false, ticks=false)
            # hidedecorations!(ax3, label=false, ticklabels=false, ticks=false)
            # hidedecorations!(ax4, label=false, ticklabels=false, ticks=false)
            # maxAe = []
            # scatterlines!(ax1, resDf.k, resDf.xbi)
            # scatterlines!(ax2, resDf.k, resDf.dbi)
            # scatterlines!(ax3, resDf.k, resDf.si)
            # scatterlines!(ax4, resDf.k, resDf.chi)

            # Label(f[:, 0], text="Score", rotation=pi/2)
            # Label(f[end+1, :], text="k")

            # CairoMakie.save(fn, f)

            return clustData, Matrix(permutedims(embedding)), cor .* exp.(-1 .* abs.(log.(abs.(grad))))
        end
        export SpectralClustering
    ##### SpectralClustering #####

    ##### draw network #####
        """
        argument
        df: dataframe of correlation matrix (return value of CCM())
        clust: cluster number of each gene (one of return value of SpectralClustering())
        pos: gene position for drawing network (one of return value of SpectralClustering())
        il: module number list which you want to draw

        return
        nw: undirected graph
        pos: node position
        cnctdf: converted correlation matrix
        clust: cluster number of each gene in network
        score: node scores
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
        fn: figure save name
        nw: network graph (one of return value of SetNetwork())
        pos: node position (one of return value of SetNetwork())
        cnctdf: converted correlation matrix (one of return value of SetNetwork())
        clust: cluster number of each gene in network (one of return value of SetNetwork())
        k: number of clusters
        node_scores: weight of node
        node_labels: label of node
        node_scaler: multiple for node diameter adjustment (Default: 100)
        edge_mode: mode of edges to be drawn
                :ALL -> All edges are drawn (Default)
                :N -> Only draw edges with negative values.
                :P -> Only draw edges with positive values.
        edge_threshold: Threshold value of edges to be drawn. (Default: 0.0)
        edge_scaler: multiple for edge thickness adjustment (Default: 5)
        x_size, y_size: Size of the drawing area (Default: 100, 100)
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
                color_list = range(LCHuv(65,100,15), stop=LCHuv(65,100,375), length=k)

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
        function DrawHeatMap(Ax::Axis, X::Vector{Int64}, Y::Vector{Int64}, Mat::Matrix{Float64}, Cnt::Matrix{Int64}, FontSize::Int64)
            heatmap!(Ax, X, Y, Mat, colormap=:bwr, colorrange=(-1.0,1.0))
            if (FontSize > 0)
                for i in X
                    for j in Y
                        txt = string(round(Mat[i,j],RoundNearestTiesAway,digits=3)) * " (n=" * string(Cnt[j]) * ")"
                        text!(Ax, position=(X[i],Y[j]), txt, fontsize=5)
                    end
                end
            end
            return nothing
        end
        """
        argument
        df1: dataframe of gene expression
        df2: dataframe of Phenomenon
        clust: cluster number of each gene (one of return value of SpectralClustering())
        fn: fig save name
        cor_mode: all gene average (:A_AVG) / positive correlation gene average (:P_AVG) / negative correlation gene average (:N_AVG) / All three types are drawn (:ALL; default)
        FontSize: set when drawing the average value of the correlation coefficient on the heat map
        """
        function CorPhenMod(df1::DataFrame, df2::DataFrame, clust::Vector{Int64}, fn::String; cor_mode::Symbol=:ALL, FontSize::Int64=0)
            kmax = length(unique(clust))

            CorList = []
            for i in 1:ncol(df2)
                push!(CorList, [[] for i=1:kmax])
            end
            for k in 1:kmax
                for i in 1:ncol(df2)
                    buf = df1[clust .== k, :]
		            Result_each = zeros(nrow(buf))
                    for j in 1:nrow(buf)
                        Result_each[j] = cor(Array(buf[j,:]), df2[:,i])
                    end
		            CorList[i][k] = Result_each
                end
            end

            x = collect(1:ncol(df2))
            y = collect(kmax:-1:1)

            f = Figure(size=(ncol(df2)*400+200, kmax*30))
            ax = []
            for i in 1:length(CorList)
                push!(ax, Axis(f[1, i], xgridvisible=false, ygridvisible=false, xticksvisible=false, yticksvisible=false, xticks=collect(-1.0:0.5:1.0)))
                if (i==1)
                    ax[1].yticks = (y,"module ".*string.(collect(1:kmax)))
                else
                    ax[i].yticklabelsvisible = false
                    linkyaxes!(ax[1], ax[i])
                end
                hidespines!(ax[i])
                ax[i].title = names(df2)[i]
                for k in 1:length(CorList[1])
                    density!(ax[i], convert.(Float64,CorList[i][k]), bins=100, scale_to=0.1, offset=(kmax-k+1), color=:x, colormap=(:bwr,0.4), colorrange=(-1.0,1.0), strokewidth=1, strokecolor=:black)
                end
            end
            save(fn, f)
        end
        export CorPhenMod
    ##### Correlation of Phenomenon and Modules #####
end

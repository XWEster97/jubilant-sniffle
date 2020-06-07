ENV["GKS_ENCODING"] = "utf-8"
include("../src/MRLSBCG.jl")
include("../src/DiagonalCG.jl")
using DataFrames,Plots,JLD

function test_run_MRLSBCG(m,n,d)
    @assert (m,d)<(n,n) " m,d<n"
    A = rand(Uniform(1,10000),(m,n))#Matrix{Float64}(I, 100, 100)
    X = zeros(n,d)
    B = rand(Uniform(-1000,1000),(m,d))
# b=reshape(Vector{Float64}(1:2000),100,20)
    pro=@timed  MRLSBCG!(X, A, B,log=true,
            verbose=true,maxiter=400,
            tol=sqrt(eps(Float64(1000))))#*d))
    return vcat(pro[1],collect(pro[2:5]))
end

function test_run_DCG(m,n,d)
    @assert (m,d)<(n,n) " m,d<n"
    A = rand(Uniform(1,10000),(m,n))#Matrix{Float64}(I, 100, 100)
    X = zeros(n,d)
    B = rand(Uniform(-1000,1000),(m,d))
# b=reshape(Vector{Float64}(1:2000),100,20)
    pro=@timed  DiagonalCG!(X, A, B,log=true,
            verbose=true,maxiter=2000,
            tol=sqrt(eps(Float64(1000))))#*d))
    return vcat(pro[1],collect(pro[2:5]))
end
# c,d=test_run(10,20,5)
# typeof(c)
function test_reserve!(proformance,pro)
    # println(size(proformance),size([pro]))
    push!(proformance, pro)
    # println(size(proformance),length([pro]))
end

function test_loop(t)

    MRLSBCGPerformance = DataFrame(history=[],size=[],
                t=[], bytes=[], gctime=[], memallocs=[])
    DiagonalPerformance=DataFrame(history=[],size=[],
                t=[], bytes=[], gctime=[], memallocs=[])
    for b in [i^2 for i=5:5:t]
        println("****----------------------d=",b,"--------MRLSBCG-----****")
        pro=test_run_MRLSBCG(6*b+10,8*b+1,Int(ceil(0.3b)))
        test_reserve!(MRLSBCGPerformance,pro)
        b=4*Int(sqrt(b))
        println("****-------DiagonalCG------d=",b,"-------------------****")
        pro=test_run_DCG(6*b+10,8*b+1,Int(ceil(0.3b)))
        test_reserve!(DiagonalPerformance,pro)
    end

    return MRLSBCGPerformance,DiagonalPerformance
end
# test_loop(9)
# result=test_loop(40)
b=40*6
@timev test_run_MRLSBCG(6*b+10,8*b+1,Int(ceil(0.3b)))
# println("****************************************","Saving")
# save("../data/mytestgit.jld", "MRLSBCG", result[1],"Diagonal",result[2])
# println("****************************************","Finished")

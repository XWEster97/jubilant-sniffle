ENV["GKS_ENCODING"] = "utf-8"
include("../src/MRLSBCG.jl")
include("../src/DiagonalCG.jl")

using DataFrames,Plots,JLD

function test_run_MRLSBCG(m,n,d;
    A = rand(Uniform(1,10000),(m,n)),
    X = zeros(n,d),
    B = rand(Uniform(-1000,1000),(m,d))
    )
    @assert (m,d)<(n,n) " m,d<n"
    pro=@timed  MRLSBCG!(X, A, B,log=true,
            verbose=true,maxiter=400,
            tol=sqrt(eps(Float64(1000)))*1000)#*d))
    return vcat(pro[1],collect(pro[2:5]))
end

function test_run_DCG(m,n,d;
    A = rand(Uniform(1,10000),(m,n)),
    X = zeros(n,d),
    B = rand(Uniform(-1000,1000),(m,d))
    )
    @assert (m,d)<(n,n) " m,d<n"
# b=reshape(Vector{Float64}(1:2000),100,20)
    pro=@timed  DiagonalCG!(X, A, B,log=true,
            verbose=true,maxiter=2000,
            tol=sqrt(eps(Float64(1000)))*1000)
    return vcat(pro[1],collect(pro[2:5]))
end

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
            # println(pro)
    end
    # println(typeof(pro))
    return MRLSBCGPerformance,DiagonalPerformance
end
test_loop(9)   ### PRECOMPILE
result=test_loop(45) ##performance test
result[1]
result[2]
println("****************************************","Performance Saving")
save("data/performancerun1.jld", "MRLSBCG", result[1],"Diagonal",result[2])
println("****************************************","Performance Finished")

function ConvergenceTest(m,n,d)
    ConMRLSBCGPerformance = DataFrame(history=[],size=[],
                t=[], bytes=[], gctime=[], memallocs=[])
    ConDiagonalPerformance=DataFrame(history=[],size=[],
                t=[], bytes=[], gctime=[], memallocs=[])
    A = rand(Uniform(1,10000),(m,n))#Matrix{Float64}(I, 100, 100)
    X = zeros(n,d)
    B = rand(Uniform(-1000,1000),(m,d))
    pro2=  test_run_DCG(m,n,d)
    test_reserve!(ConDiagonalPerformance,pro2)
    pro2=0#free memory
    pro1= test_run_MRLSBCG(m,n,d)
    test_reserve!(ConMRLSBCGPerformance,pro1)
    return ConMRLSBCGPerformance,ConDiagonalPerformance
end
# a=ConvergenceTest(25,50,5)

conresult=ConvergenceTest(1300,2000,180)

println("****************************************","Recording Convergence Test")
save("data/maccorun.jld", "MRLSBCG", conresult[1],"Diagonal",conresult[2])
println("****************************************","Finished All")

####################test
# b=45*6
# @timev test_run_MRLSBCG(6*b+10,8*b+1,Int(ceil(0.5b)))

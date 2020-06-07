import SparseArrays, IterativeSolvers, Base.iterate
using  Random, Distributions, LinearAlgebra, ProgressMeter, Profile,BenchmarkTools
SA,IS=SparseArrays,IterativeSolvers
using Printf


mutable struct DiagonalCGIterable{matT, solT, vecT, numT <: Real}
    A::matT
    x::solT
    r::vecT
    c::vecT
    u::vecT
    reltol::numT
    residual::numT
    prev_residual::numT
    maxiter::Int
    mv_products::Int
end

mutable struct MRLSBCGIterable{matT, solT,numT <: Real}
    A::matT
    X::solT
    G::matT
    S::matT
    P::matT
    Q::matT
    reltol::numT
    residual::numT
    prev_G::matT
    maxiter::Int
    mv_products::Int
end

@inline converged(it::Union{MRLSBCGIterable,DiagonalCGIterable}) = it.residual ≤ it.reltol

@inline start(it::Union{MRLSBCGIterable,DiagonalCGIterable}) = 0

@inline done(it::Union{MRLSBCGIterable,DiagonalCGIterable}, iteration::Int) = iteration ≥ it.maxiter || converged(it)

function generate_SPD(n)
    A=rand(DiscreteUniform(-100,100),(n,n))
    D=rand(Uniform(1,100),n)
    return A*Diagonal(D)*A'
end

breshape(B)=reshape(B,size(B,1)*size(B,2),1)

function CGblockdiag(d,X::SA.SparseMatrixCSC)
    X=fill(X,d)
    print(size(X))
    num = length(X)
    mX = Int[ size(x, 1) for x in X ]
    nX = Int[ size(x, 2) for x in X ]
    m = sum(mX)
    n = sum(nX)
    X[1].nzval
    Tv = promote_type(map(x->eltype(x.nzval), X)...)
    Ti = isempty(X) ? Int : promote_type(map(x->eltype(x.rowval), X)...)

    colptr = Vector{Ti}(undef, n+1)
    nnzX = Int[ SA.nnz(x) for x in X ]
    nnz_res = sum(nnzX)
    rowval = Vector{Ti}(undef, nnz_res)
    nzval = Vector{Tv}(undef, nnz_res)

    let nnz_sofar = 0 ;nX_sofar = 0; mX_sofar = 0
        for i = 1 : num
            colptr[(1 : nX[i] + 1) .+ nX_sofar] = X[i].colptr .+ nnz_sofar
            rowval[(1 : nnzX[i]) .+ nnz_sofar] = X[i].rowval .+ mX_sofar
            nzval[(1 : nnzX[i]) .+ nnz_sofar] = X[i].nzval
            nnz_sofar += nnzX[i]
            nX_sofar += nX[i]
            mX_sofar += mX[i]
        end
        colptr[n+1] = nnz_sofar + 1
    end
    SA.SparseMatrixCSC(m, n, colptr, rowval, nzval)
end

#coding=utf-8
Include("head.jl")

struct MRLSBCGStateVariables{T,Tx<:AbstractArray{T}}
    P::Tx
    G::Tx
    Q::Tx
end

function iterate(it::MRLSBCGIterable, iteration::Int=start(it))
    mul!(it.Q,it.A,it.P)
    #αQQSS
    # SS=it.S'*it.S
    # print(it.Q)
    SS=it.S'*it.S
    α= try
        inv(it.Q'*it.Q)*SS
    catch y
        pinv(it.Q'*it.Q)*SS
    end
    # α=inv(it.Q'*it.Q)*it.S'*it.S
    #XXPα
    # α=inv(it.Q'*it.Q)*SS
    it.X +=it.P*α
    it.G -=it.Q*α
    it.residual = norm(it.G )
    if done(it, iteration)
        return nothing
    end
    SS_P=SS
    mul!(it.S,it.A',it.G)
    SS=it.S'*it.S
    β=inv(SS_P)*SS
    it.P=it.S+it.P*β
    # print(norm(it.G),"\n")
    norm(it.G), iteration + 1
end

function MRLSBCG_iterator!(X, A, B, Pl = IS.Identity();
    tol = sqrt(eps(real(eltype(B))))*size(b,2)*size(b,1),
    maxiter::Int = size(A, 2)*size(B,2)+100,
    statevars::MRLSBCGStateVariables = MRLSBCGStateVariables(zero(X), similar(X), similar(X)),
    initially_zero::Bool = false)
    P = statevars.P
    G = statevars.G
    Q = statevars.Q
    G = B
    # copyto!(G, B)
    mv_products = 0
    S = A'*G
    Q = A* S
    P = S
    SS=S'*S
    residual = norm(G)
    reltol = residual * tol
    @printf("%s\t%1.2e\n","reltol=",reltol)
    if isa(Pl, IS.Identity)
        return MRLSBCGIterable(A, X, G, S, P, Q,
            reltol, residual, zero(G),
            maxiter, mv_products
        )
    else
        return PCGIterable(Pl, A, x, r, c, u,
            reltol, residual, one(eltype(x)),
            maxiter, mv_products
        )
    end
end

function MRLSBCG!(X, A, B;
    tol = sqrt(eps(real(eltype(B))))*size(b,2)*size(b,1),
    maxiter::Int = size(A, 2)*size(B,2),
    log::Bool = false,
    statevars::MRLSBCGStateVariables = MRLSBCGStateVariables(zero(X), similar(X), similar(X)),
    verbose::Bool = false,
    Pl = IS.Identity(),
    kwargs...
)

    history = IS.ConvergenceHistory(partial = !log)
    history[:tol] = tol
    log && IS.reserve!(history, :resnorm, maxiter + 1)

    # Actually perform MRLSBCG
    iterable = MRLSBCG_iterator!(X, A, B, Pl; tol = tol, maxiter = maxiter, statevars = statevars, kwargs...)
    if log
        history.mvps = iterable.mv_products
    end
    for (iteration, item) = enumerate(iterable)
        if log
            IS.nextiter!(history, mvps = 1)
            push!(history, :resnorm, iterable.residual)
        end
        verbose && @printf("%3d\t%1.2e\n", iteration, iterable.residual)
    end

    verbose && println()
    log && IS.setconv(history, converged(iterable))
    log && IS.shrink!(history)

    log ? [history, [size(A,1),size(A,2),size(B,2)]] : iterable.X

end
# A= rand(Uniform(1,10000),(1000,1800))#Matrix{Float64}(I, 100, 100)
# X=zeros(1800,36)
# B=rand(Uniform(-1000,1000),(1000,36))
# # b=reshape(Vector{Float64}(1:2000),100,20)
# @timev ch= cg!(X, A, B,verbose=true,maxiter=1000,log=true,tol=sqrt(eps(real(mean(B)))))
# eps(mean(B))
#
# # typeof(ch)
# ch[1]
# # df=DataFrame()
# # push!(df,[1 2 3])
# norm(B)*sqrt(eps(real(mean(B))))
# ch[1]

#coding=utf-8
include("./head.jl")

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

function iterate(it::DiagonalCGIterable, iteration::Int=start(it))
    if done(it, iteration) return nothing end
    # u := r + βu (almost an axpy)
    β = it.residual^2 / it.prev_residual^2
    it.u .= it.r .+ β .* it.u
    # c = A * u
    mul!(it.c, it.A, it.u)
    α = it.residual^2 / dot(it.u, it.c)

    # Improve solution and residual
    it.x .+= α .* it.u
    it.r .-= α .* it.c

    it.prev_residual = it.residual
    it.residual = norm(it.r)

    # Return the residual at item and iteration number as state
    it.residual, iteration + 1
end

function Diagonal_cg_iterator!(x, A, b, Pl = IS.Identity();
    tol = sqrt(eps(real(eltype(b)))),
    maxiter::Int = size(A, 2),
    statevars::IS.CGStateVariables = IS.CGStateVariables(zero(x), similar(x), similar(x)),
    initially_zero::Bool = false
)
    u = statevars.u
    r = statevars.r
    c = statevars.c
    u .= zero(eltype(x))
    copyto!(r, b)

	mv_products = 0
	c = similar(x)
	residual = norm(b)
	reltol = residual * tol
    DiagonalCGIterable(A, x, r, c, u,
    reltol, residual, one(residual),
    maxiter, mv_products)
end



function DiagonalCG!(x, A, B;
    tol = sqrt(eps(real(eltype(B)))),
    maxiter::Int = size(A, 2),
    log::Bool = false,
    statevars::IS.CGStateVariables = IS.CGStateVariables(zero(x), similar(x), similar(x)),
    verbose::Bool = false,
    Pl = IS.Identity(),
    kwargs...)
    b=breshape(A'*B)
    A=CGblockdiag(size(B,2),SA.sparse(A'*A))
    x=breshape(x)
    statevars=IS.CGStateVariables(zero(x), similar(x), similar(x))
    history = IS.ConvergenceHistory(partial = !log)
    history[:tol] = tol
    log && IS.reserve!(history, :resnorm, maxiter + 1)

    iterable = Diagonal_cg_iterator!(x, A, b, Pl; tol = tol, maxiter = maxiter, statevars = statevars, kwargs...)
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
	SB=size(B,2)
    log ? [history, [size(A,1)/SB,size(A,2)/SB,SB]] : iterable.x
end

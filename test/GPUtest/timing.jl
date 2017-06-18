function mv_dd(A::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64}; maxiter=20)
    d_A = CudaSparseMatrixCSR(A)
    d_x = CudaArray(x)
    d_y = CudaArray(Float64, size(A,1))
    for count = 1 : maxiter
        d_y = CUSPARSE.mv('N',d_A, d_x,'0')
    end
    y = to_host(d_y)
    return y
end


function hybmv(elty,m,n,arr)
    A = convert(SparseMatrixCSC{elty,Int},sprand(m,n,prob))
    x = rand(elty,n)
    y = rand(elty,m)
    alpha = rand(elty)
    beta  = rand(elty)
    tic()
    d_A = CudaSparseMatrixCSR(A)
    d_A = CUSPARSE.switch2hyb(d_A)
    d_x = CudaArray(x)
    d_y = CudaArray(y)
    for i in 1:20
        d_y = CUSPARSE.hybmv!('N',alpha,d_A,d_x,beta,d_y,'O')
    end
    h_y = to_host(d_y)
    ti = toq()
    return vcat(arr,[ti])
end

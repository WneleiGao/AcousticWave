using ParSpMatVec

n = 1000000;
sp= 1e-7;
A = sprand(n,n,sp);
x = randn(n);
y = zeros(n);

@time AmulB!(y,A,x);
@time y1 = A * x;

N = 20; tcost = zeros(N);
for i = 1 : N
    tic()
    AmulB!(y,A,x);
    tcost[i] = toq();
end
mean(tcost)

N = 20; tcost1 = zeros(N);
for i = 1 : N
    tic()
    A_mul_B!(y1, A, x)
    tcost1[i] = toq();
end
mean(tcost1) / mean(tcost)
norm(y-y1)

# ===================test complex number ====================
n = 1000000;
sp= 1e-7;
A = sprand(n,n,sp);
x = randn(n) + im*randn(n);
y = zeros(Complex128, n);
y1 = zeros(Complex128, n);

@time AmulB!(y,A,x);
@time A_mul_B!(y1,A,x);

N = 20; tcost = zeros(N);
for i = 1 : N
    tic()
    AmulB!(y,A,x);
    tcost[i] = toq();
end

N = 20; tcost1 = zeros(N);
for i = 1 : N
    tic()
    A_mul_B!(y1,A,x)
    tcost1[i] = toq();
end
mean(tcost1) / mean(tcost)
norm(y1-y)

# ===================test complex matrix times complex vector====================
n = 1000000;
sp= 1e-7;
A = sprand(Complex128, n,n,sp)
x = randn(n) + im*randn(n);
y = zeros(Complex128, n);
y1 = zeros(Complex128, n);

@time AmulB!(y,A,x);
@time A_mul_B!(y1,A,x);

N = 20; tcost = zeros(N);
for i = 1 : N
    tic()
    AmulB!(y,A,x);
    tcost[i] = toq();
end

N = 20; tcost1 = zeros(N);
for i = 1 : N
    tic()
    A_mul_B!(y1,A,x)
    tcost1[i] = toq();
end
mean(tcost1) / mean(tcost)
norm(y1-y)


# ============================adjoint========================================

n = 1000000;
sp= 1e-7;
A = sprand(n,n,sp);
x = randn(n);
y = zeros(n);

@time AcmulB!(y,A,x);
@time At_mul_B!(y1, A, x);

N = 20; tcost = zeros(N);
K = 10;
tcost = zeros(N, K)
for n = 1 : K
    for i = 1 : N
        tic()
        AcmulB!(y,A,x, nthreads=n);
        tcost[i, n] = toq();
    end
end

N = 20; tcost1 = zeros(N);
for i = 1 : N
    tic()
    At_mul_B!(y1, A, x)
    tcost1[i] = toq();
end
mean(tcost1) ./ mean(tcost,1)
norm(y-y1)

# ===================test complex number ====================
n = 1000000;
sp= 1e-7;
A = sprand(n,n,sp);
x = randn(n)+im*randn(n);
y = zeros(Complex128, n);

@time AcmulB!(y,A,x);
@time At_mul_B!(y1, A, x);

N = 20; tcost = zeros(N);
K = 10;
tcost = zeros(N, K)
for n = 1 : K
    for i = 1 : N
        tic()
        AcmulB!(y,A,x, nthreads=n);
        tcost[i, n] = toq();
    end
end

N = 20; tcost1 = zeros(N);
for i = 1 : N
    tic()
    At_mul_B!(y1, A, x)
    tcost1[i] = toq();
end
mean(tcost1) ./ mean(tcost,1)
norm(y-y1)
# ===================test complex matrix times complex vector====================
n = 1000000;
sp= 1e-7;
A = sprand(Complex128, n,n,sp);
x = randn(n) + im*randn(n);
y = zeros(Complex128, n);
y1 = zeros(Complex128, n);

@time AcmulB!(y,A,x);
@time Ac_mul_B!(y1, A, x);

N = 20; tcost = zeros(N);
K = 10;
tcost = zeros(N, K)
for n = 1 : K
    for i = 1 : N
        tic()
        AcmulB!(y,A,x, nthreads=n);
        tcost[i, n] = toq();
    end
end

N = 20; tcost1 = zeros(N);
for i = 1 : N
    tic()
    Ac_mul_B!(y1, A, x)
    tcost1[i] = toq();
end
mean(tcost1) ./ mean(tcost,1)
norm(y-y1)

type FidMtx
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     dz :: Float64
     dx :: Float64
     dt :: Float64
     MvzBvz :: Array{Float64,1}
     MvzBp  :: SparseMatrixCSC{Float64,Int64}
     MvxBvx :: Array{Float64,1}
     MvxBp  :: SparseMatrixCSC{Float64,Int64}
     MpzBpz :: Array{Float64,1}
     MpzBvz :: SparseMatrixCSC{Float64,Int64}
     MpxBpx :: Array{Float64,1}
     MpxBvx :: SparseMatrixCSC{Float64,Int64}
end

type gpuFdMtx
     MvzBvz :: CUSPARSE.CudaSparseMatrixHYB{Float32}
     MvzBp  :: CUSPARSE.CudaSparseMatrixHYB{Float32}
     MvxBvx :: CUSPARSE.CudaSparseMatrixHYB{Float32}
     MvxBp  :: CUSPARSE.CudaSparseMatrixHYB{Float32}
     MpzBpz :: CUSPARSE.CudaSparseMatrixHYB{Float32}
     MpzBvz :: CUSPARSE.CudaSparseMatrixHYB{Float32}
     MpxBpx :: CUSPARSE.CudaSparseMatrixHYB{Float32}
     MpxBvx :: CUSPARSE.CudaSparseMatrixHYB{Float32}
end

function SendFidMtx2GPU(fm::FidMtx)
 gm=gpuFdMtx(CUSPARSE.switch2hyb(CudaSparseMatrixCSR(spdiagm(convert(Vector{Float32}, fm.MvzBvz)))),
             CUSPARSE.switch2hyb(CudaSparseMatrixCSR(convert(SparseMatrixCSC{Float32,Int64}, fm.MvzBp))),
             CUSPARSE.switch2hyb(CudaSparseMatrixCSR(spdiagm(convert(Vector{Float32}, fm.MvxBvx)))),
             CUSPARSE.switch2hyb(CudaSparseMatrixCSR(convert(SparseMatrixCSC{Float32,Int64}, fm.MvxBp))),
             CUSPARSE.switch2hyb(CudaSparseMatrixCSR(spdiagm(convert(Vector{Float32}, fm.MpzBpz)))),
             CUSPARSE.switch2hyb(CudaSparseMatrixCSR(convert(SparseMatrixCSC{Float32,Int64}, fm.MpzBvz))),
             CUSPARSE.switch2hyb(CudaSparseMatrixCSR(spdiagm(convert(Vector{Float32}, fm.MpxBpx)))),
             CUSPARSE.switch2hyb(CudaSparseMatrixCSR(convert(SparseMatrixCSC{Float32,Int64}, fm.MpxBvx))))
    return gm
end

function DispStable!(vmax::Float64, vmin::Float64, f0::Float64, dz::Float64, dx::Float64, dt::Float64)
    h = vmin / 5 / (f0*2)
    d = minimum([dz, dx])
    tt = 6*d / (7*sqrt(3)*vmax)
    if dx >= h || dz >= h || dt >= tt
       println("maximum acceptable grid size: $h")
       println("maximum acceptable time step: $tt")
       error("unstable or frequency dispersion")
    end
    return nothing
end

function InitFidMtx(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dz::Float64, dx::Float64, dt::Float64, vmin::Float64, vmax::Float64, v::Array{Float64,2}, f0::Float64)
    DispStable!(vmax, vmin, f0, dz, dx, dt)
    Cpml = InitCpml(nz, nx, ext, vmax, dz, dx, iflag)
    # Cpml = InitCpml(nz, nx, ext, vmax, rho, dz, dx, iflag)
    (MvzBvz, MvzBp)  = Mvz(Cpml.vz,    dz, dt, ext, iflag)
    (MvxBvx, MvxBp)  = Mvx(Cpml.vx,    dx, dt, ext, iflag)
    (MpzBpz, MpzBvz) = Mpz(Cpml.pz, v, dz, dt, ext, iflag)
    (MpxBpx, MpxBvx) = Mpx(Cpml.px, v, dx, dt, ext, iflag)
    fidMtx = FidMtx(nz, nx, ext, iflag, dz, dx, dt, MvzBvz, MvzBp, MvxBvx, MvxBp, MpzBpz, MpzBvz, MpxBpx, MpxBvx)
    return fidMtx
end




function vp2lambda(vp::Array{Float64,2}, rho::Array{Float64,2})
    (m , n ) = size(vp)
    (m1, n1) = size(rho)
    if m1 != m || n1 != n
       error("size not match")
    end
    lambda = Array(typeof(vp[1]), m, n)
    for j = 1 : n
        for i = 1 : m
            lambda[i,j] = rho[i,j] * vp[i,j]^2
        end
    end
    return lambda
end

function lambda2vp(lambda::Array{Float64, 2}, rho::Array{Float64, 2})
    (m, n) = size(lambda)
    (m1, n1) = size(rho)
    if m1 != m || n1 != n
       error("size not match")
    end
    vp = Array(typeof(lambda[1]), m, n)
    for j = 1 : n
        for i = 1 : m
            vp[i,j] = sqrt(lambda[i,j] / rho[i,j])
        end
    end
    return vp
end

function Mvz(vz_pml::SparseMatrixCSC{Float64,Int64}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
    (m,n) = size(vz_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dz; c2 = a2/dz;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + vz_pml[iz,ix]/2
            b = 1/dt - vz_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = 1 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(m, m)
    tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
    for iz = 2: m-2
        tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
        tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
    end
    tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
    tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
    MvzBvz = C3
    MvzBp  = spdiagm(denum) * kron(speye(n), tmp)
    return MvzBvz, MvzBp
end

function Mvx(vx_pml::SparseMatrixCSC{Float64,Int64}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
    (m,n) = size(vx_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dx; c2 = a2/dx;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + vx_pml[iz,ix]/2
            b = 1/dt - vx_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = 1 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(n, n)
    tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
    for ix = 2: n-2
        tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
        tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
    end
    tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
    tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
    MvxBvx = C3
    MvxBp  = spdiagm(denum) * kron(tmp, speye(m))
    return MvxBvx, MvxBp
end

function Mpz(pz_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
    v = modExpand(v, ext, iflag)
    (m, n) = size(pz_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dz; c2 = a2/dz;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + pz_pml[iz,ix]/2
            b = 1/dt - pz_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(m, m)
    tmp[1,1] =  c1; tmp[1,2] = c2;
    tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
    for iz = 3: m-1
        tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
        tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
    end
    tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
    MpzBpz = C3
    MpzBvz = spdiagm(denum) * kron(speye(n), tmp)
    return MpzBpz, MpzBvz
end

function Mpx(px_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
    v = modExpand(v, ext, iflag)
    (m, n) = size(px_pml)
    a1 = 9/8  ; a2 = -1/24;
    c1 = a1/dx; c2 = a2/dx;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + px_pml[iz,ix]/2
            b = 1/dt - px_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(n, n)
    tmp[1,1] =  c1; tmp[1,2] = c2;
    tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
    for ix = 3: n-1
        tmp[ix,ix-1] = -c1; tmp[ix, ix  ] = c1;
        tmp[ix,ix-2] = -c2; tmp[ix, ix+1] = c2;
    end
    tmp[n,n-1] = -c1; tmp[n,n] = c1; tmp[n,n-2] = -c2;
    MpxBpx = C3
    MpxBvx = spdiagm(denum) * kron(tmp, speye(m))
    return MpxBpx, MpxBvx
end

# the finite difference matrix for back ward propagation
# function InitFidMtx_back(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dz::Float64, dx::Float64, dt::Float64, vmax::Float64, vmin::Float64, f0::Float64, v::Array{Float64,2})
#     DispStable!(vmax, vmin, f0, dz, dx, dt)
#     Cpml = InitCpml(nz, nx, ext, vmax, dz, dx, iflag)
#     # Cpml = InitCpml(nz, nx, ext, vmax, rho, dz, dx, iflag)
#     (MvzBvz, MvzBp)  = Mvz_back(Cpml.vz,    dz, dt, ext, iflag)
#     (MvxBvx, MvxBp)  = Mvx_back(Cpml.vx,    dx, dt, ext, iflag)
#     (MpzBpz, MpzBvz) = Mpz_back(Cpml.pz, v, dz, dt, ext, iflag)
#     (MpxBpx, MpxBvx) = Mpx_back(Cpml.px, v, dx, dt, ext, iflag)
#     fidMtx = FidMtx(nz, nx, ext, iflag, dz, dx, dt, MvzBvz, MvzBp, MvxBvx, MvxBp, MpzBpz, MpzBvz, MpxBpx, MpxBvx)
#     return fidMtx
# end

# function Mpx(px_pml::SparseMatrixCSC{Float64,Int64}, lambda::Array{Float64,2}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     lambda = modExpand(lambda, ext, iflag)
#     (m, n) = size(lambda)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dx; c2 = a2/dx;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt + px_pml[iz,ix]/2
#             b = 1/dt - px_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = lambda[iz,ix] / a
#             C3[(ix-1)*m+iz]= b / a
#         end
#     end
#     tmp = spzeros(n, n)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#     for ix = 3: n-1
#         tmp[ix,ix-1] = -c1; tmp[ix, ix  ] = c1;
#         tmp[ix,ix-2] = -c2; tmp[ix, ix+1] = c2;
#     end
#     tmp[n,n-1] = -c1; tmp[n,n] = c1; tmp[n,n-2] = -c2;
#     MpxBpx = spdiagm(C3)
#     MpxBvx = spdiagm(denum) * kron(tmp, speye(m))
#     return  MpxBpx, MpxBvx
# end

# function Mpz(pz_pml::SparseMatrixCSC{Float64,Int64}, lambda::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     lambda = modExpand(lambda, ext, iflag)
#     (m, n) = size(lambda)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dz; c2 = a2/dz;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt + pz_pml[iz,ix]/2
#             b = 1/dt - pz_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = lambda[iz,ix] / a
#             C3[(ix-1)*m+iz]= b / a
#         end
#     end
#     tmp = spzeros(m, m)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#     for iz = 3: m-1
#         tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
#         tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
#     end
#     tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
#     MpzBpz = spdiagm(C3)
#     MpzBvz = spdiagm(denum) * kron(speye(n), tmp)
#     return MpzBpz, MpzBvz
# end


# function Mpz_back(pz_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     v = modExpand(v, ext, iflag)
#     (m, n) = size(pz_pml)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dz; c2 = a2/dz;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt - pz_pml[iz,ix]/2
#             b = 1/dt + pz_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
#             C3[(ix-1)*m+iz]= b / a
#         end
#     end
#     tmp = spzeros(m, m)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#     for iz = 3: m-1
#         tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
#         tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
#     end
#     tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
#     MpzBpz = spdiagm(C3)
#     MpzBvz = -spdiagm(denum) * kron(speye(n), tmp)
#     return MpzBpz, MpzBvz
# end

# function Mvx(vx_pml::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     rho = modExpand(rho, ext, iflag)
#     (m,n) = size(rho)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dx; c2 = a2/dx;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             if ix < n
#                a = (rho[iz,ix+1]+rho[iz,ix]) / (2*dt)
#             else
#                a = rho[iz,ix] / dt
#             end
#             b = vx_pml[iz,ix] / 2
#             denum[(ix-1)*m+iz] = 1 / (a+b)
#             C3[(ix-1)*m+iz]= (a-b) / (a+b)
#         end
#     end
#     tmp = spzeros(n, n)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#     for ix = 2: n-2
#         tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
#         tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
#     end
#     tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
#     tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
#     MvxBvx = spdiagm(C3)
#     MvxBp = spdiagm(denum) * kron(tmp, speye(m))
#     return MvxBvx, MvxBp
# end

# function Mvx_back(vx_pml::SparseMatrixCSC{Float64,Int64}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     (m,n) = size(vx_pml)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dx; c2 = a2/dx;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt - vx_pml[iz,ix]/2
#             b = 1/dt + vx_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = 1 / a
#             C3[(ix-1)*m+iz]= b / a
#         end
#     end
#     tmp = spzeros(n, n)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#     for ix = 2: n-2
#         tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
#         tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
#     end
#     tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
#     tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
#     MvxBvx = spdiagm(C3)
#     MvxBp  = -spdiagm(denum) * kron(tmp, speye(m))
#     return MvxBvx, MvxBp
# end

# function Mvz_back(vz_pml::SparseMatrixCSC{Float64,Int64}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     (m,n) = size(vz_pml)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dz; c2 = a2/dz;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt - vz_pml[iz,ix]/2
#             b = 1/dt + vz_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = 1 / a
#             C3[(ix-1)*m+iz]= b / a
#         end
#     end
#     tmp = spzeros(m, m)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#     for iz = 2: m-2
#         tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
#         tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
#     end
#     tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
#     tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
#     MvzBvz = spdiagm(C3)
#     MvzBp  = -spdiagm(denum) * kron(speye(n), tmp)
#     return MvzBvz, MvzBp
# end

# function Mvz(vz_pml::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     rho = modExpand(rho, ext, iflag)
#     (m,n) = size(rho)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dz; c2 = a2/dz;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             if iz < m
#                a = (rho[iz+1,ix]+rho[iz,ix]) / (2*dt)
#             else
#                a = rho[iz,ix] / dt
#             end
#             b = vz_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = 1 / (a+b)
#             C3[(ix-1)*m+iz]= (a-b) / (a+b)
#         end
#     end
#     tmp = spzeros(m, m)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#     for iz = 2: m-2
#         tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
#         tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
#     end
#     tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
#     tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
#     MvzBvz = spdiagm(C3)
#     MvzBp  = spdiagm(denum) * kron(speye(n), tmp)
#     return MvzBvz, MvzBp
# end

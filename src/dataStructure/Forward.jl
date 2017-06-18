function sup1!(c::Array{Float64,1}, a::Array{Float64,1}, b::Array{Float64,1})
    @inbounds for i = 1 : length(a)
        c[i] = a[i]*b[i]
    end
    return nothing
end

function sup2!(c::Array{Float64,1}, a::Array{Float64,1}, b::Array{Float64,1})
    @inbounds for i = 1 : length(a)
        c[i] = a[i]+b[i]
    end
    return nothing
end

function sup3!(a::Array{Float64,1}, b::Array{Float64,1})
    @inbounds for i = 1 : length(a)
        a[i] = a[i]+b[i]
    end
    return nothing
end

function OneStepForward!(spt2::SnapShot, spt1::SnapShot, fidMtx::FidMtx)
    spt2.it  = spt1.it + 1
    spt2.vz = fidMtx.MvzBvz .* spt1.vz + fidMtx.MvzBp  * (spt1.pz+spt1.px)
    spt2.vx = fidMtx.MvxBvx .* spt1.vx + fidMtx.MvxBp  * (spt1.pz+spt1.px)
    spt2.pz = fidMtx.MpzBpz .* spt1.pz + fidMtx.MpzBvz *  spt2.vz
    spt2.px = fidMtx.MpxBpx .* spt1.px + fidMtx.MpxBvx *  spt2.vx
    return nothing
end

function OneStepForwardGPU!(spt2::gpuSpt, spt1::gpuSpt, GM::gpuFdMtx,
                            tmp::CUDArt.CudaArray{Float32,1})
    n = length(spt1.px);
    spt2.it = spt1.it + 1;
    CUSPARSE.mv!('N', 1.f0, GM.MvzBvz, spt1.vz, 0.0f0, spt2.vz, 'O');
    CUBLAS.blascopy!(n, spt1.px, 1, tmp, 1);
    CUBLAS.axpy!(n, 1.f0, spt1.pz, 1, tmp, 1);
    CUSPARSE.mv!('N', 1.f0, GM.MvzBp, tmp, 1.0f0, spt2.vz, 'O');

    CUSPARSE.mv!('N', 1.f0, GM.MvxBvx, spt1.vx, 0.0f0, spt2.vx, 'O');
    CUSPARSE.mv!('N', 1.f0, GM.MvxBp, tmp, 1.0f0, spt2.vx, 'O');

    CUSPARSE.mv!('N', 1.f0, GM.MpzBpz, spt1.pz, 0.0f0, spt2.pz, 'O');
    CUSPARSE.mv!('N', 1.f0, GM.MpzBvz, spt2.vz, 1.0f0, spt2.pz, 'O');

    CUSPARSE.mv!('N', 1.f0, GM.MpxBpx, spt1.px, 0.0f0, spt2.px, 'O');
    CUSPARSE.mv!('N', 1.f0, GM.MpxBvx, spt2.vx, 1.0f0, spt2.px, 'O');

    return nothing
end

function OneStepForward!(spt2::SnapShot, spt1::SnapShot, fidMtx::FidMtx, tmp::Array{Float64,1}, tmp1::Array{Float64,1})
    spt2.it  = spt1.it + 1
    sup1!(spt2.vz, fidMtx.MvzBvz, spt1.vz)
    sup2!(tmp, spt1.pz, spt1.px)
    AmulB!(tmp1, fidMtx.MvzBp, tmp)
    sup3!(spt2.vz, tmp1)

    sup1!(spt2.vx, fidMtx.MvxBvx, spt1.vx)
    AmulB!(tmp1, fidMtx.MvxBp, tmp)
    sup3!(spt2.vx, tmp1)

    sup1!(spt2.pz, fidMtx.MpzBpz, spt1.pz)
    AmulB!(tmp1, fidMtx.MpzBvz, spt2.vz)
    sup3!(spt2.pz, tmp1)

    sup1!(spt2.pz, fidMtx.MpzBpz, spt1.pz)
    AmulB!(tmp1, fidMtx.MpzBvz, spt2.vz)
    sup3!(spt2.pz, tmp1)

    sup1!(spt2.px, fidMtx.MpxBpx, spt1.px)
    AmulB!(tmp1, fidMtx.MpxBvx, spt2.vx)
    sup3!(spt2.px, tmp1)
    return nothing
end

# output one shot, inject windowed source cube
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, wsc::WSrcCube, fidMtx::FidMtx; tmax=1.0)
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt)+1
    shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddWsc2spt!(spt1, wsc)
    spt2shot!(shot, spt1)
    tmp = zeros(length(spt1.vz))
    tmp1= zeros(tmp);
    for it = 2 : wsc.nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        AddWsc2spt!(spt2, wsc)
        copySnapShot!(spt1, spt2)
        spt2shot!(shot, spt1)
    end
    for it = wsc.nt+1 : nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        copySnapShot!(spt1, spt2)
        spt2shot!(shot, spt1)
    end
    return shot
end

# output shot or shot3, inject single source
# function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, src::Source, fidMtx::FidMtx; tmax=1.0, otype="p")
#     nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
#     nt = round(Int64, tmax/dt)+1
#     isz= src.isz; isx= src.isx;
#     if otype == "p"
#        shot = InitShot(isz, isx, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     elseif otype == "vp"
#        shot = InitShot3(isz, isx, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
#     end
#     tl = src.ot; tu = tl + (src.nt-1)*dt;
#     spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
#     AddSource!(spt1, src)
#     tmp = zeros(length(spt1.vz))
#     tmp1= zeros(tmp);
#     if otype == "p"
#        spt2shot!(shot, spt1)
#     elseif otype == "vp"
#        spt2shot3!(shot, spt1)
#     end
#     for it = 2: nt
#         OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
#         if tl <= (it-1)*dt <= tu
#            AddSource!(spt2, src)
#         end
#         copySnapShot!(spt1, spt2)
#         if otype == "p"
#            spt2shot!(shot, spt1)
#         elseif otype == "vp"
#            spt2shot3!(shot, spt1)
#         end
#     end
#     return shot
# end

function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, src::Source, fidMtx::FidMtx; tmax=1.0, otype="p")
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt)+1
    isz= src.isz; isx= src.isx;
    if otype == "p"
       shot = InitShot(isz, isx, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    elseif otype == "vp"
       shot = InitShot3(isz, isx, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    end
    tl = src.ot; tu = tl + (src.nt-1)*dt;
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddSource!(spt1, src)
    if otype == "p"
       spt2shot!(shot, spt1)
    elseif otype == "vp"
       spt2shot3!(shot, spt1)
    end
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx)
        if tl <= (it-1)*dt <= tu
           AddSource!(spt2, src)
        end
        copySnapShot!(spt1, spt2)
        if otype == "p"
           spt2shot!(shot, spt1)
        elseif otype == "vp"
           spt2shot3!(shot, spt1)
        end
    end
    return shot, spt1
end


function MultiStepForwardGPU(irz::Vector{Int64}, irx::Vector{Int64}, src::Source, gm::gpuFdMtx; tmax=1.0)
    nt = floor(Int64, tmax/src.dt) + 1
    nz = src.nz; nx = src.nx; ext = src.ext;
    iflag = src.iflag; dt = src.dt;
    if iflag == 1
       Nz = nz + 2*ext
       zupper = ext
    else src.iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    N  = Nz * Nx
    s1 = InitGpuSpt(nz, nx, ext, iflag, dt, 1)
    s2 = InitGpuSpt(nz, nx, ext, iflag, dt, 2)
    p  = convert(Vector{Float32}, src.p)
    idx = (src.isx+ext-1)*Nz + src.isz + zupper
    scatterSrc!(s1, p, N, idx)
    tmp = CudaArray(zeros(Float32, N))
    for i = 2 : length(p)
        OneStepForwardGPU!(s2, s1, gm, tmp)
        scatterSrc!(s2, p, N, idx)
        copySptGPU!(s1, s2)
    end
    for i = length(p)+1 : nt
        OneStepForwardGPU!(s2, s1, gm, tmp)
        copySptGPU!(s1, s2)
    end
    return s1
end

# output shot or shot3, inject multi-sources
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, srcs::Array{Source,1}, fidMtx::FidMtx; tmax=1.0, otype="p")
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt)+1
    (tl, tu) = SrcRange(srcs)
    if otype == "p"
       shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    elseif otype == "vp"
       shot = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    end
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddMultiSources!(spt1, srcs)
    tmp = zeros(length(spt1.vz))
    tmp1= zeros(tmp);
    if otype == "p"
       spt2shot!(shot, spt1)
    elseif otype == "vp"
       spt2shot3!(shot, spt1)
    end
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if tl <= (it-1)*dt <= tu
           AddMultiSources!(spt2, srcs)
        end
        copySnapShot!(spt1, spt2)
        if otype == "p"
           spt2shot!(shot, spt1)
        elseif otype == "vp"
           spt2shot3!(shot, spt1)
        end
    end
    return shot
end

# save snapshot // full wave field // stress field to slow memory, inject single source
function MultiStepForward(path::String , src::Source, fidMtx::FidMtx; tmax=1.0, wtype="p")
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt)+1
    tl = src.ot; tu = tl + (src.nt-1)*dt;
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddSource!(spt1, src)
    tmp = zeros(length(spt1.vz))
    tmp1= zeros(tmp);
    if wtype == "spt"
       fid = writeSnapShot(path, spt1)
    elseif wtype == "vp"
       fid = writeWfd(path, spt1)
    elseif wtype == "p"
       fid = writeStress(path, spt1)
    end
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if tl <= (it-1)*dt <= tu
           AddSource!(spt2, src)
        end
        copySnapShot!(spt1, spt2)
        if wtype == "spt"
           writeSnapShot(fid, spt1)
        elseif wtype == "vp"
           writeWfd(fid, spt1)
        elseif wtype == "p"
           writeStress(fid, spt1)
        end
    end
    close(fid)
    return nothing
end

# save full wave field or stress field to slow memory, inject multi-sources
function MultiStepForward(path::String , srcs::Array{Source,1}, fidMtx::FidMtx; tmax=1.0, wtype="p")
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt)+1
    (tl, tu) = SrcRange(srcs)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddMultiSources!(spt1, srcs)
    tmp = zeros(length(spt1.vz))
    tmp1= zeros(tmp);
    if wtype == "spt"
       fid = writeSnapShot(path, spt1)
    elseif wtype == "vp"
       fid = writeWfd(path, spt1)
    elseif wtype == "p"
       fid = writeStress(path, spt1)
    end
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if tl <= (it-1)*dt <= tu
           AddMultiSources!(spt2, srcs)
        end
        copySnapShot!(spt1, spt2)
        if wtype == "spt"
           writeSnapShot(fid, spt1)
        elseif wtype == "vp"
           writeWfd(fid, spt1)
        elseif wtype == "p"
           writeStress(fid, spt1)
        end
    end
    close(fid)
    return nothing
end

# return shot or shot3, inject stress which read from slow memory
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, path::String , fidMtx::FidMtx; tmax=1.0, otype="p")
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt)+1
    (nz1, nx1, ext1, iflag1, dt1, nt1) = InfoStress(path)
    if otype == "p"
       shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    elseif otype == "vp"
       shot = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    end
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    addStress2Spt!(spt1, path)
    tmp = zeros(length(spt1.vz))
    tmp1= zeros(tmp)
    if otype == "p"
       spt2shot!(shot, spt1)
    elseif otype == "vp"
       spt2shot3!(shot, spt1)
    end
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if it <= nt1
           addStress2Spt!(spt2, path)
        end
        copySnapShot!(spt1, spt2)
        if otype == "p"
           spt2shot!(shot, spt1)
        elseif otype == "vp"
           spt2shot3!(shot, spt1)
        end
    end
    return shot
end

# return shot or shot3, do born approximation
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, I::Array{Float64,2}, path::String , fidMtx::FidMtx; tmax=1.0, otype="p")
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt) + 1
    (nz1, nx1, ext1, iflag1, dt1, nt1) = InfoStress(path)
    if otype == "p"
       shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    elseif otype == "vp"
       shot = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    end
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    addReflection2Spt!(spt1, I, path)
    tmp = zeros(length(spt1.vz))
    tmp1= zeros(tmp);
    if otype == "p"
       spt2shot!(shot, spt1)
    elseif otype == "vp"
       spt2shot3!(shot, spt1)
    end
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if it <= nt1
           addReflection2Spt!(spt2, I, path)
        end
        copySnapShot!(spt1, spt2)
        if otype == "p"
           spt2shot!(shot, spt1)
        elseif otype == "vp"
           spt2shot3!(shot, spt1)
        end
    end
    return shot
end

#  output shot or shot3, input wavelet and the distribution of source(used for source location)
function MultiStepForward(irz::Array{Int64,1}, irx::Array{Int64,1}, w::Array{Float64,1}, dis::Array{Float64,2}, fidMtx::FidMtx; tmax=1.0, otype="p")
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt) + 1
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nx +   ext
    end
    Nx  = nx + 2*ext
    nst = length(w)
    if otype == "p"
       shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    elseif otype == "vp"
       shot = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    end
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddMultiSources!(spt1, w, dis)
    tmp = zeros(length(spt1.vz))
    tmp1= zeros(tmp);
    if otype == "p"
       spt2shot!(shot, spt1)
    elseif otype == "vp"
       spt2shot3!(shot, spt1)
    end
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1)
        if it <= nst
           AddMultiSources!(spt2, w, dis)
        end
        copySnapShot!(spt1, spt2)
        if otype == "p"
           spt2shot!(shot, spt1)
        elseif otype == "vp"
           spt2shot3!(shot, spt1)
        end
    end
    return shot
end

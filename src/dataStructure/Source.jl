type Source
     isz :: Int64
     isx :: Int64
     nz  :: Int64
     nx  :: Int64
     ext :: Int64
     iflag :: Int64
     ot  :: Float64
     dt  :: Float64
     nt  :: Int64
     p   :: Array{Float64,1}
end

function Ricker(f0::Float64, dt::Float64)
	  nw = 2.2/f0/dt
	  nw = 2*floor(Int,nw/2)+1
	  nc = floor(Int,nw/2)
	  w  = zeros(nw)
	  k  = collect(1:nw)
    k  = vec(k)
	  alpha = (nc-k+1)*f0*dt*pi
	  beta = alpha.^2
	  w = (1.-beta.*2).*exp(-beta)
	  return w
end

function randAmp(ns::Int64)
    amp = zeros(ns)
    for is = 1 : ns
        if rand() <=0.5
           amp[is] = -1.0
        else
           amp[is] =  1.0
        end
    end
    return amp
end

function randST(ns::Int64, dt::Float64, s::Float64)
    ot = s .* randn(ns)
    ts = minimum(ot)
    ot = (ot + abs(ts))/2.0
    indt = round(Int64, ot/dt)
    ot[:] = indt*dt
    return ot
end

function InitSource(isz::Int64, isx::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, ot::Float64, dt::Float64, f0::Float64, amp::Float64)
    p = Ricker(f0, dt)
    p = amp * p
    nt= length(p)
    src = Source(isz, isx, nz, nx, ext, iflag, ot, dt, nt, p)
    return src
end

function InitMultiSources(isz::Array{Int64,1}, isx::Array{Int64,1}, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, ot::Array{Float64,1}, dt::Float64, f0::Float64, amp::Array{Float64,1})
    ns = length(isz)
    srcs = Array(Source, ns)
    for is = 1 : ns
        srcs[is] = InitSource(isz[is], isx[is], nz, nx, ext, iflag, ot[is], dt, f0, amp[is])
    end
    return srcs
end

function scatterSrc!(gs::gpuSpt, p::Vector{Float32}, N::Int64, idx::Int64)
    it = gs.it
    s = sparsevec([idx], p[it], N)
    ds= CudaSparseVector(s)
    CUSPARSE.sctr!(ds, gs.pz, 'O')
    return nothing
end

function AddSource!(spt::SnapShot, src::Source)
    dt = spt.dt
    nz = spt.nz ; nx    = spt.nx   ;
    ext= spt.ext; iflag = spt.iflag;
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
    end
    Nx  = nx + 2*ext
    it  = spt.it
    if src.ot <= (it-1)*dt <= src.ot+(src.nt-1)*dt
       indt = it - round(Int64, src.ot/dt)
       indz = src.isz + zupper
       indx = src.isx + ext
       pos = (indx-1) * Nz + indz
       spt.pz[pos] = spt.pz[pos] + (src.p[indt])/2
       spt.px[pos] = spt.px[pos] + (src.p[indt])/2
    end
    return nothing
end

# give multi-sources
function AddMultiSources!(spt::SnapShot, srcs::Array{Source,1})
    nz=spt.nz; nx=spt.nx; ext=spt.ext; iflag=spt.iflag; dt=spt.dt
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
    end
    Nx  = nx + 2*ext
    it  = spt.it
    ns  = length(srcs)
    for is = 1: ns
        if srcs[is].ot <= (it-1)*dt <= srcs[is].ot+(srcs[is].nt-1)*dt
           indt = it - round(Int64, srcs[is].ot/dt)
           indz = srcs[is].isz + zupper
           indx = srcs[is].isx + ext
           pos = (indx-1) * Nz + indz
           spt.pz[pos] = spt.pz[pos] + srcs[is].p[indt]/2
           spt.px[pos] = spt.px[pos] + srcs[is].p[indt]/2
        end
    end
    return nothing
end

#give distribution and wavelet of sources
function AddMultiSources!(spt::SnapShot, w::Array{Float64,1}, dis::Array{Float64,2})
    it = spt.it
    if it > length(w)
       error("it is out the range of wavelet")
    end
    nz = spt.nz ; nx = spt.nx;
    ext= spt.ext; iflag = spt.iflag;
    if length(dis) != nz*nx
       error("size of dis dismatch with the size of SnapShot")
    end
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
       Nx  = nx + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
       Nx  = nx + 2*ext
    end
    tmp = zeros(Nz, Nx);
    tmp[zupper+1: zupper+nz, ext+1:ext+nx] = 1/2 * w[it] * dis
    tmp = vec(tmp)
    spt.pz[:] = spt.pz[:] + tmp[:]
    spt.px[:] = spt.px[:] + tmp[:]
    return nothing
end

function SrcRange(srcs::Array{Source,1})
    tmax = 0.0
    tmin = 0.0
    ns   = length(srcs)
    dt   = srcs[1].dt
    for is = 1: ns
        tmp = srcs[is].ot
        tmp1= tmp + (srcs[is].nt-1)*dt
        if tmp < tmin
           tmin = tmp
        end
        if tmp1 > tmax
           tmax = tmp1
        end
    end
    return tmin, tmax
end

function srcs2spt(srcs::Array{Source,1}, it::Int64)
    nz = srcs[1].nz ; nx = srcs[1].nx;
    ext= srcs[1].ext; iflag=srcs[1].iflag;
    dt = srcs[1].dt
    spt = InitSnapShot(nz, nx, ext, iflag, dt, it)
    AddMultiSources!(spt, srcs)
    return spt
end

function WriteSrcs2SnapShots(path::String , srcs::Array{Source,1})
    nz = srcs[1].nz
    (tl, tu) = SrcRange(srcs)
    nt  = round(Int64, (tu-tl)/dt) + 1
    fid = open(path, "w")
    write(fid, nz, nx, ext, iflag, dt)
    for it = 1 : nt
        it1 = round(Int64, (tl+(it-1)*dt)/dt+1)
        spt = srcs2spt(srcs, it1)
        writeSnapShot(fid, spt)
    end
    close(fid)
    return nothing
end

function writeSrcs2Stress(path::String , srcs::Array{Source,1})
    ns = length(srcs);
    nz = srcs[1].nz; nx = srcs[1].nx; ext = srcs[1].ext; iflag = srcs[1].iflag; dt = srcs[1].dt
    (tl, tu) = SrcRange(srcs)
    nt  = round(Int64, (tu-tl)/dt) + 1
    fid = open(path, "w")
    write(fid, nz, nx, ext, iflag, dt)
    for it = 1 : nt
        p   = zeros(nz, nx)
        for is = 1 : ns
            if srcs[is].ot <= tl+(it-1)*dt <= srcs[is].ot + dt*(srcs[is].nt-1)
               indt = it - round(Int64, (srcs[is].ot-tl)/dt)
               isz = srcs[is].isz; isx = srcs[is].isx;
               p[isz,isx] = p[isz,isx] + srcs[is].p[indt]
            end
        end
        write(fid, vec(p))
    end
    close(fid)
    return nothing
end

function multiStress2wlet(path::String , iz::Int64, ix::Int64)
    fid = open(path, "r")
    (nz, nx, ext, iflag, dt, nt) = InfoStress(path)
    w = zeros(nt)
    for it = 1 : nt
        sts= readStress(path, it)
        w[it] = sts.p[iz,ix]
    end
    return w
end

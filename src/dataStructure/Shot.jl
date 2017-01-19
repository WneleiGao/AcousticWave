type Shot
     nr  :: Int64
     isz :: Int64
     isx :: Int64
     nz  :: Int64
     nx  :: Int64
     ext :: Int64
     iflag::Int64
     irz :: Array{Int64,1}
     irx :: Array{Int64,1}
     ot  :: Float64
     dt  :: Float64
     nt  :: Int64
     p   :: Array{Float64,2}
end

type Shot3
     nr  :: Int64
     isz :: Int64
     isx :: Int64
     nz  :: Int64
     nx  :: Int64
     ext :: Int64
     iflag::Int64
     irz :: Array{Int64,1}
     irx :: Array{Int64,1}
     ot  :: Float64
     dt  :: Float64
     nt  :: Int64
     vz  :: Array{Float64,2}
     vx  :: Array{Float64,2}
     p   :: Array{Float64,2}
end

function InitShot(isz::Int64, isx::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, irz::Array{Int64,1}, irx::Array{Int64,1}, ot::Float64, dt::Float64, nt::Int64)
    nr = length(irz)
    s  = Shot(nr, isz, isx, nz, nx, ext, iflag, irz, irx, ot, dt, nt, zeros(nt,nr))
    return s
end

function InitShot3(isz::Int64, isx::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, irz::Array{Int64,1}, irx::Array{Int64,1}, ot::Float64, dt::Float64, nt::Int64)
    nr = length(irz)
    s3 = Shot3(nr, isz, isx, nz, nx, ext, iflag, irz, irx, ot, dt, nt, zeros(nt,nr), zeros(nt,nr), zeros(nt,nr))
    return s3
end

function writeShot!(path::String , shot::Shot)
    fid = open(path, "w")
    nr  = length(shot.irz)
    write(fid, Int32(nr),      Int32(shot.isz), Int32(shot.isx))
    write(fid, Int32(shot.nz), Int32(shot.nx) , Int32(shot.ext), Int32(shot.iflag))
    write(fid, convert(Array{Int32},shot.irz) , convert(Array{Int32},shot.irx))
    write(fid, Float32(shot.ot), Float32(shot.dt), Int32(shot.nt))
    write(fid, convert(Array{Float32},vec(shot.p)))
    close(fid)
    return nothing
end

function writeShot!(path::String , shot::Shot3)
    fid = open(path, "w")
    nr  = length(shot.irz)
    write(fid, Int32(nr),      Int32(shot.isz), Int32(shot.isx))
    write(fid, Int32(shot.nz), Int32(shot.nx) , Int32(shot.ext), Int32(shot.iflag))
    write(fid, convert(Array{Int32},shot.irz) , convert(Array{Int32},shot.irx))
    write(fid, Float32(shot.ot), Float32(shot.dt), Int32(shot.nt))
    write(fid, convert(Array{Float32},vec(shot.vz)), convert(Array{Float32},vec(shot.vx)), convert(Array{Float32},vec(shot.p)))
    close(fid)
    return nothing
end

function readShot(path::String )
    fid = open(path, "r")
    nr  = Int64(read(fid,Int32)); isz = Int64(read(fid,Int32)); isx = Int64(read(fid,Int32));
    nz  = Int64(read(fid,Int32)); nx  = Int64(read(fid,Int32)); ext = Int64(read(fid,Int32)); iflag = Int64(read(fid,Int32));
    irz = convert(Array{Int64},read(fid,Int32,nr)); irx = convert(Array{Float64},read(fid,Int32,nr));
    ot  = Float64(read(fid,Float32)); dt = Float64(read(fid,Float32)); nt = Int64(read(fid,Int32));
    p   = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    close(fid)
    shot = Shot(nr, isz, isx, nz, nx, ext, iflag, irz, irx, ot, dt, nt, p)
    return shot
end

function readShot3(path::String )
    fid = open(path, "r")
    nr  = Int64(read(fid,Int32)); isz = Int64(read(fid,Int32)); isx = Int64(read(fid,Int32));
    nz  = Int64(read(fid,Int32)); nx  = Int64(read(fid,Int32)); ext = Int64(read(fid,Int32)); iflag = Int64(read(fid,Int32));
    irz = convert(Array{Int64},read(fid,Int32,nr)); irx = convert(Array{Float64},read(fid,Int32,nr));
    ot  = Float64(read(fid,Float32)); dt = Float64(read(fid,Float32)); nt = Int64(read(fid,Int32));
    vz  = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    vx  = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    p   = convert(Array{Float64}, reshape(read(fid,Float32,nt*nr),nt,nr))
    close(fid)
    shot = Shot3(nr, isz, isx, nz, nx, ext, iflag, irz, irx, ot, dt, nt, vz, vx, p)
    return shot
end

function spt2shot!(shot::Shot, spt::SnapShot)
    it = spt.it;
    nz = shot.nz; nx=shot.nx; ext=shot.ext; iflag=shot.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       zupper = ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    nr = shot.nr
    for ir = 1 : nr
        iz = shot.irz[ir] + zupper
        ix = shot.irx[ir] + ext
        ind= (ix-1)*Nz + iz
        shot.p[it, ir] = spt.pz[ind] + spt.px[ind]
    end
    return nothing
end

function spt2shot3!(shot::Shot3, spt::SnapShot)
    it = spt.it;
    nz=shot.nz; nx=shot.nx; ext=shot.ext; iflag=shot.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       zupper = ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    nr = shot.nr
    for ir = 1 : nr
        iz = shot.irz[ir] + zupper
        ix = shot.irx[ir] + ext
        ind= (ix-1)*Nz + iz
        shot.vz[it, ir]= spt.vz[ind]
        shot.vx[it, ir]= spt.vx[ind]
        shot.p[it, ir] = spt.pz[ind] + spt.px[ind]
    end
    return nothing
end

function AddShot2SnapShot!(spt::SnapShot, shot::Shot)
    it = spt.it;
    nz = shot.nz; nx=shot.nx; ext=shot.ext; iflag=shot.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       zupper = ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    nr = shot.nr
    for ir = 1 : nr
        iz = shot.irz[ir] + zupper
        ix = shot.irx[ir] + ext
        ind= (ix-1)*Nz + iz
        spt.pz[ind] = spt.pz[ind] + shot.p[it,ir]
        spt.px[ind] = spt.px[ind] + shot.p[it,ir]
    end
end

function AddShot32SnapShot!(spt::SnapShot, shot::Shot3)
    it = spt.it;
    nz = shot.nz; nx=shot.nx; ext=shot.ext; iflag=shot.iflag;
    if iflag == 1
       Nz = nz + 2*ext
       zupper = ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    nr = shot.nr
    for ir = 1 : nr
        iz = shot.irz[ir] + zupper
        ix = shot.irx[ir] + ext
        ind= (ix-1)*Nz + iz
        spt.vz[ind] = spt.vz[ind] + shot.vz[it, ir]
        spt.vx[ind] = spt.vx[ind] + shot.vx[it, ir]
        spt.pz[ind] = spt.pz[ind] + shot.p[it,  ir]
        spt.px[ind] = spt.px[ind] + shot.p[it,  ir]
    end
end

function copyShot(s::Shot)
    s1 = Shot(s.nr, s.isz, s.isx, s.nz, s.nx, s.ext, s.iflag, s.irz, s.irx, s.ot, s.dt, s.nt, s.p)
    return s1
end

function copyShot3(s::Shot3)
    s1 = Shot3(s.nr, s.isz, s.isx, s.nz, s.nx, s.ext, s.iflag, s.irz, s.irx, s.ot, s.dt, s.nt, s.vz, s.vx, s.p)
    return s1
end

function wfd2Shot(irz::Array{Int64,1}, irx::Array{Int64,1}, path::String )
    (nz, nx, ext, iflag, dt, nt) = InfoWfd(path)
    nr  = length(irz)
    shot= InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    fid = open(path, "r")
    for it = 1 : nt
        position = sizeof(Float64)*5 + (it-1)*nz*nx*sizeof(Float64)*3 + nz*nx*sizeof(Float64)*2
        seek(fid, position)
        p = reshape(read(fid, Float64, nz*nx), nz, nx)
        for ir = 1 : nr
            iz = shot.irz[ir]
            ix = shot.irx[ir]
            shot.p[it,ir] = p[iz, ix]
        end
    end
    close(fid)
    return shot
end

function wfd2Shot3(irz::Array{Int64,1}, irx::Array{Int64,1}, path::String )
    (nz, nx, ext, iflag, dt, nt) = InfoWfd(path)
    nr  = length(irz)
    shot= InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    fid = open(path, "r")
    for it = 1 : nt
        position = sizeof(Float64)*5 + (it-1)*nz*nx*sizeof(Float64)*3
        seek(fid, position)
        vz = reshape(read(fid, Float64, nz*nx), nz, nx)
        vx = reshape(read(fid, Float64, nz*nx), nz, nx)
        p  = reshape(read(fid, Float64, nz*nx), nz, nx)
        for ir = 1 : nr
            iz = shot.irz[ir]
            ix = shot.irx[ir]
            shot.vz[it,ir] = vz[iz, ix]
            shot.vx[it,ir] = vx[iz, ix]
            shot.p[it,ir]  = p[iz, ix]
        end
    end
    close(fid)
    return shot
end

function stress2Shot(irz::Array{Int64,1}, irx::Array{Int64,1}, path::String )
    (nz, nx, ext, iflag, dt, nt) = InfoStress(path)
    nr  = length(irz)
    shot= InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    fid = open(path, "r")
    for it = 1 : nt
        position = sizeof(Float64)*5 + (it-1)*nz*nx*sizeof(Float64)
        seek(fid, position)
        p = reshape(read(fid, Float64, nz*nx), nz, nx)
        for ir = 1 : nr
            iz = shot.irz[ir]
            ix = shot.irx[ir]
            shot.p[it,ir] = p[iz, ix]
        end
    end
    close(fid)
    return shot
end

function minusShot(shot1::Shot, shot2::Shot)
    if shot1.irz != shot2.irz || shot1.irx != shot2.irx || shot1.nt != shot2.nt
       error("size dismatch")
    end
    shot = InitShot(shot1.isz, shot1.isx, shot1.nz, shot1.nx, shot1.ext, shot1.iflag, shot1.irz, shot1.irx, 0.0, shot1.dt, shot1.nt)
    shot.d[:] = shot1.d - shot2.d
    return shot
end

function minusShot3(shot1::Shot3, shot2::Shot3)
    if shot1.irz != shot2.irz || shot1.irx != shot2.irx || shot1.nt != shot2.nt
       error("size dismatch")
    end
    shot = InitShot(shot1.isz, shot1.isx, shot1.nz, shot1.nx, shot1.ext, shot1.iflag, shot1.irz, shot1.irx, shot1.dt, shot1.nt)
    shot.vz[:] = shot1.vz - shot2.vz
    shot.vx[:] = shot1.vx - shot2.vx
    shot.p[:] = shot1.p - shot2.p
    return shot
end

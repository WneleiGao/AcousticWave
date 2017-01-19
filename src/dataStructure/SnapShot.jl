type SnapShot
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     dt :: Float64
     it :: Int64
     vz :: Array{Float64,1}
     vx :: Array{Float64,1}
     pz :: Array{Float64,1}
     px :: Array{Float64,1}
end

type Wfd
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     dt :: Float64
     it :: Int64
     vz :: Array{Float64,2}
     vx :: Array{Float64,2}
     p  :: Array{Float64,2}
end

type Stress
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     dt :: Float64
     it :: Int64
     p  :: Array{Float64,2}
end

function InitSnapShot(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    if iflag == 1
       Nz = nz + 2*ext
       Nx = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       Nx = nx + 2*ext
    end
    spt = SnapShot(nz, nx, ext, iflag, dt, it, zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx))
    return spt
end

function InitWfd(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    spt = Wfd(nz, nx, ext, iflag, dt, it, zeros(nz,nx), zeros(nz,nx), zeros(nz,nx))
    return spt
end

function InitStress(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    spt = Stress(nz, nx, ext, iflag, dt, it, zeros(nz,nx))
    return spt
end

function writeSnapShot(path::String , spt::SnapShot)
    fid = open(path, "w")
    write(fid, spt.nz ); write(fid, spt.nx   );
    write(fid, spt.ext); write(fid, spt.iflag);
    write(fid, spt.dt );
    write(fid, spt.vz ); write(fid, spt.vx   );
    write(fid, spt.pz ); write(fid, spt.px   );
    flush(fid)
    return fid
end

function writeSnapShot(fid::IOStream, spt::SnapShot)
    write(fid, spt.vz); write(fid, spt.vx);
    write(fid, spt.pz); write(fid, spt.px);
    flush(fid)
    return nothing
end

function writeWfd(path::String , spt::SnapShot)
    fid = open(path, "w")
    nz = spt.nz; nx = spt.nx; ext= spt.ext; iflag = spt.iflag; dt = spt.dt;
    write(fid, nz, nx, ext, iflag, dt);
    if spt.iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif spt.iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    tmp = reshape(spt.vz, Nz, Nx);
    tmp1= tmp[zupper+1:zupper+nz, ext+1:ext+nx]
    write(fid, vec(tmp1))
    tmp = reshape(spt.vx, Nz, Nx);
    tmp1= tmp[zupper+1:zupper+nz, ext+1:ext+nx]
    write(fid, vec(tmp1))
    tmp = spt.pz + spt.px; tmp = reshape(tmp, Nz, Nx);
    tmp1= tmp[zupper+1:zupper+nz, ext+1:ext+nx]
    write(fid, vec(tmp1))
    flush(fid)
    return fid
end

function writeWfd(fid::IOStream, spt::SnapShot)
    nz = spt.nz; nx = spt.nx; ext= spt.ext; iflag = spt.iflag;
    if iflag == 1
       zupper = spt.ext
       Nz = spt.nz + 2*spt.ext
    elseif iflag == 2
       Nz = spt.nz +   spt.ext
       zupper = 0
    end
    Nx = spt.nx + 2*spt.ext
    tmp = reshape(spt.vz, Nz, Nx);
    tmp1= tmp[zupper+1:zupper+nz, ext+1:ext+nx]
    write(fid, vec(tmp1))
    tmp = reshape(spt.vx, Nz, Nx);
    tmp1= tmp[zupper+1:zupper+nz, ext+1:ext+nx]
    write(fid, vec(tmp1))
    tmp = spt.pz + spt.px; tmp = reshape(tmp, Nz, Nx);
    tmp1= tmp[zupper+1:zupper+nz, ext+1:ext+nx]
    write(fid, vec(tmp1))
    flush(fid)
    return nothing
end

function writeStress(path::String , spt::SnapShot)
    fid = open(path, "w")
    nz = spt.nz; nx = spt.nx; ext= spt.ext; iflag = spt.iflag; dt = spt.dt;
    write(fid, nz, nx, ext, iflag, dt)
    if spt.iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif spt.iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    tmp = spt.pz + spt.px; tmp = reshape(tmp, Nz, Nx);
    tmp1= tmp[zupper+1:zupper+nz, ext+1:ext+nx]
    write(fid, vec(tmp1))
    flush(fid)
    return fid
end

function writeStress(fid::IOStream, spt::SnapShot)
    nz = spt.nz; nx = spt.nx; ext= spt.ext; iflag = spt.iflag;
    if spt.iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif spt.iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    tmp = spt.pz + spt.px; tmp = reshape(tmp, Nz, Nx);
    tmp1= tmp[zupper+1:zupper+nz, ext+1:ext+nx]
    write(fid, vec(tmp1))
    flush(fid)
    return nothing
end

function readSnapShot(path::String , it::Int64)
    fid = open(path, "r") ;
    nz  = read(fid, Int64); nx    = read(fid, Int64);
    ext = read(fid, Int64); iflag = read(fid, Int64);
    dt  = read(fid, Float64);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    position = sizeof(Int64)*5 + (it-1)*sizeof(Float64)*Nz*Nx*4
    seek(fid, position)
    vz = read(fid, Float64, Nz*Nx)
    vx = read(fid, Float64, Nz*Nx)
    pz = read(fid, Float64, Nz*Nx)
    px = read(fid, Float64, Nz*Nx)
    spt = SnapShot(nz, nx, ext, iflag, dt, it, vz, vx, pz, px)
    close(fid)
    return spt
end

function readWfd(path::String , it::Int64)
    fid = open(path, "r") ;
    nz  = read(fid, Int64); nx = read(fid, Int64);
    ext = read(fid, Int64); iflag = read(fid, Int64);
    dt  = read(fid, Float64);
    position = sizeof(Int64)*5 + (it-1)*sizeof(Float64)*nz*nx*3
    seek(fid, position)
    vz = reshape(read(fid, Float64, nz*nx), nz, nx)
    vx = reshape(read(fid, Float64, nz*nx), nz, nx)
    p  = reshape(read(fid, Float64, nz*nx), nz, nx)
    wfd= Wfd(nz, nx, ext, iflag, dt, it, vz, vx, p)
    close(fid)
    return wfd
end

function readStress(path::String , it::Int64)
    fid = open(path, "r") ;
    nz  = read(fid, Int64); nx = read(fid, Int64);
    ext = read(fid, Int64); iflag = read(fid, Int64);
    dt  = read(fid, Float64);
    position = sizeof(Int64)*5 + (it-1)*sizeof(Float64)*nz*nx
    seek(fid, position)
    p   = reshape(read(fid, Float64, nz*nx), nz, nx)
    sts = Stress(nz, nx, ext, iflag, dt, it, p)
    close(fid)
    return sts
end

function InfoSnapShot(path::String ; print_flag = false)
    fid = open(path, "r")
    nz  = read(fid, Int64); nx  = read(fid, Int64);
    ext = read(fid, Int64); iflag=read(fid, Int64);
    dt =read(fid, Float64);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    nt = round(Int64, (filesize(fid) - sizeof(Int64)*5) / (4*Nz*Nx*sizeof(Float64)))
    close(fid)
    if print_flag
       println("zlength: $nz, xlength: $nx")
       println("padding: $ext, surface: $iflag")
       println("dt: $dt, nt: $nt")
    end
    return nz, nx, ext, iflag, dt, nt
end

function InfoWfd(path::String ; print_flag = false)
    fid = open(path, "r")
    nz  = read(fid, Int64); nx  = read(fid, Int64);
    ext = read(fid, Int64); iflag=read(fid, Int64);
    dt =read(fid, Float64);
    nt = round(Int64, (filesize(fid) - sizeof(Int64)*5) / (3*nz*nx*sizeof(Float64)))
    close(fid)
    if print_flag
       println("zlength: $nz, xlength: $nx")
       println("padding: $ext, surface: $iflag")
       println("dt: $dt, nt: $nt")
    end
    return nz, nx, ext, iflag, dt, nt
end

function InfoStress(path::String ; print_flag = false)
    fid = open(path, "r")
    nz  = read(fid, Int64); nx  = read(fid, Int64);
    ext = read(fid, Int64); iflag=read(fid, Int64);
    dt =read(fid, Float64);
    nt = round(Int64, (filesize(fid) - sizeof(Int64)*5) / (nz*nx*sizeof(Float64)))
    close(fid)
    if print_flag
       println("zlength: $nz, xlength: $nx")
       println("padding: $ext, surface: $iflag")
       println("dt: $dt, nt: $nt")
    end
    return nz, nx, ext, iflag, dt, nt
end

function copySnapShot!(snapShot1::SnapShot, snapShot2::SnapShot)
    if snapShot2.nz != snapShot1.nz || snapShot2.nx != snapShot1.nx || snapShot2.ext != snapShot1.ext || snapShot2.dt != snapShot1.dt || snapShot2.iflag != snapShot1.iflag
       error("the two snapShot are different")
    end
    snapShot1.it  = snapShot2.it
    n = length(snapShot1.vz)
    @inbounds for i = 1 : n
        snapShot1.vz[i] = snapShot2.vz[i]
    end
    @inbounds for i = 1 : n
        snapShot1.vx[i] = snapShot2.vx[i]
    end
    @inbounds for i = 1 : n
        snapShot1.pz[i] = snapShot2.pz[i]
    end
    @inbounds for i = 1 : n
        snapShot1.px[i] = snapShot2.px[i]
    end
    return nothing
end

function normSnapShot(snapShot::SnapShot)
    L = norm(snapShot.vz)
    L = L + dot(vec(snapShot.vz), vec(snapShot.vz))
    L = L + dot(vec(snapShot.vx), vec(snapShot.vx))
    L = L + dot(vec(snapShot.pz), vec(snapShot.pz))
    L = L + dot(vec(snapShot.px), vec(snapShot.px))
    L = sqrt(L)
    return L
end

function minusSnapShot(snapShot1::SnapShot, snapShot2::SnapShot)
    nz = snapShot1.nz
    nx = snapShot1.nx
    ext = snapShot1.ext
    iflag = snapShot1.iflag
    dt = snapShot1.dt
    it = snapShot1.it
    if nz != snapShot2.nz || nx != snapShot2.nx || ext != snapShot2.ext || iflag != snapShot2.iflag
       error("size dismatch")
    end
    snapShot3 = InitSnapShot(nz, nx, ext, iflag, dt, it)
    snapShot3.vz[:] = snapShot1.vz - snapShot2.vz
    snapShot3.vx[:] = snapShot1.vx - snapShot2.vx
    snapShot3.pz[:] = snapShot1.pz - snapShot2.pz
    snapShot3.px[:] = snapShot1.px - snapShot2.px
    return snapShot3
end

function addSnapShot(snapShot1::SnapShot, snapShot2::SnapShot)
    nz = snapShot1.nz
    nx = snapShot1.nx
    ext = snapShot1.ext
    iflag = snapShot1.iflag
    dt = snapShot1.dt
    it = snapShot1.it
    if nz != snapShot2.nz || nx != snapShot2.nx || ext != snapShot2.ext || iflag != snapShot2.iflag
       error("size dismatch")
    end
    snapShot3 = InitSnapShot(nz, nx, ext, iflag, dt, it)
    snapShot3.vz[:] = snapShot1.vz + snapShot2.vz
    snapShot3.vx[:] = snapShot1.vx + snapShot2.vx
    snapShot3.pz[:] = snapShot1.pz + snapShot2.pz
    snapShot3.px[:] = snapShot1.px + snapShot2.px
    return snapShot3
end

function reverseSnapShotsOrder(path::String , path_tmp::String )
    (nz, nx, ext, iflag, dt, nt) = InfoSnapShot(path_tmp)
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    fid = open(path, "w")
    write(fid, nz, nx, ext, iflag, dt)
    fid_tmp = open(path_tmp, "r")
    for it = 1 : nt
        position = sizeof(Float64)*5 + sizeof(Float64)*Nz*Nx*4*(nt-it)
        seek(fid_tmp, position)
        d = read(fid_tmp, Float64, Nz*Nx*4);
        write(fid, d)
    end
    close(fid_tmp); close(fid)
    rm(path_tmp)
    return nothing
end

function reverseWfdOrder(path::String , path_tmp::String )
    (nz, nx, ext, iflag, dt, nt) = InfoWfd(path_tmp)
    fid = open(path, "w")
    write(fid, nz, nx, ext, iflag, dt)
    fid_tmp = open(path_tmp, "r")
    for it = 1 : nt
        position = sizeof(Float64)*5 + sizeof(Float64)*nz*nx*3*(nt-it)
        seek(fid_tmp, position)
        d = read(fid_tmp, Float64, nz*nx*3);
        write(fid, d)
    end
    close(fid_tmp); close(fid)
    rm(path_tmp)
    return nothing
end

function reverseStressOrder(path::String , path_tmp::String )
    (nz, nx, ext, iflag, dt, nt) = InfoStress(path_tmp)
    fid = open(path, "w")
    write(fid, nz, nx, ext, iflag, dt)
    fid_tmp = open(path_tmp, "r")
    for it = 1 : nt
        position = sizeof(Float64)*5 + sizeof(Float64)*nz*nx*(nt-it)
        seek(fid_tmp, position)
        d = read(fid_tmp, Float64, nz*nx);
        write(fid, d)
    end
    close(fid_tmp); close(fid)
    rm(path_tmp)
    return nothing
end

function mixNormStress(path::String )
    (nz, nx, ext, iflag, dt, nt) = InfoStress(path)
    fid = open(path, "r")
    seek(fid, sizeof(Float64)*5)
    tmp = zeros(nz*nx)
    for it = 1 : nt
        sts = read(fid, Float64, nz*nx)
        tmp = tmp + sts.^2
    end
    close(fid)
    tmp = sqrt(tmp)
    return tmp
end

function L2normStress(path::String )
    (nz, nx, ext, iflag, dt, nt) = InfoStress(path)
    fid = open(path, "r")
    seek(fid, sizeof(Float64)*5)
    tmp = 0.0
    for it = 1 : nt
        spt = read(fid, Float64, nz*nx)
        tmp = tmp + dot(spt, spt)
    end
    close(fid)
    tmp = sqrt(tmp)
    return tmp
end

function scaleStress(path::String , alpha::Float64)
    (nz, nx, ext, iflag, dt, nt) = InfoStress(path)
    fid = open(path, "r+")
    pre = sizeof(Float64) * 5
    lspt = nz*nx*sizeof(Float64)
    for it = 1 : nt
        position = pre + (it-1)*lspt
        seek(fid, position)
        spt = read(fid, Float64, nz*nx)
        spt[:] = spt * alpha
        seek(fid, position)
        write(fid, spt)
    end
    close(fid)
    return nothing
end

function addStress2Spt!(spt::SnapShot, path::String )
    it = spt.it
    nz = spt.nz;   nx = spt.nx;
    ext = spt.ext; iflag = spt.iflag;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    sts = readStress(path, it)
    p   = sts.p * 1/2
    tmp = reshape(spt.pz, Nz, Nx)
    tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + p
    spt.pz = vec(tmp)
    tmp = reshape(spt.px, Nz, Nx)
    tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + p
    spt.px = vec(tmp)
    return nothing
end

function addReflection2Spt!(spt::SnapShot, I::Array{Float64,2}, path::String )
    it = spt.it
    nz = spt.nz; nx = spt.nx;
    ext = spt.ext; iflag = spt.iflag;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    sts = readStress(path, it)
    ref = sts.p .* I * 1/2
    tmp = reshape(spt.pz, Nz, Nx)
    tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + ref
    spt.pz[:] = vec(tmp)
    tmp = reshape(spt.px, Nz, Nx)
    tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + ref
    spt.px[:] = vec(tmp)
    return nothing
end

function spt2dis(w::Array{Float64,1}, spt::SnapShot)
    nz = spt.nz ; nx = spt.nx;
    ext= spt.ext; iflag = spt.iflag;
    it = spt.it
    if length(w) < spt.it
       error("out of source time range")
    end
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    pz = reshape(spt.pz, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
    px = reshape(spt.px, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
    dis = 1/2 * w[it] * (pz+px)
    return dis
end

function spt2wlet(w::Array{Float64,1}, dis::Array{Float64,2}, spt::SnapShot)
    nz = spt.nz ; nx = spt.nx;
    ext= spt.ext; iflag = spt.iflag;
    it = spt.it
    if length(w) < it
       error("out of source time range")
    end
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    p = reshape(spt.pz, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx] + reshape(spt.px, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
    p = vec(p)
    w[it] = 1/2 * dot(vec(dis), vec(p))
    return nothing
end

function diswlet2spt(dis::Array{Float64,1}, w::Array{Float64,1}, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    if length(w) < it
       error("out of source time range")
    end
    if length(dis) != nz*nx
       error("length of dis does not much model size")
    end
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    spt = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, it)
    tmp = zeros(Nz, Nx); dis_tmp = reshape(dis, nz, nx)
    tmp[zupper+1: zupper+nz, ext+1:ext+nx] = 1/2 * w[it] * dis_tmp
    tmp = vec(tmp)
    spt.pz[:] = tmp[:]
    spt.px[:] = tmp[:]
    return spt
end

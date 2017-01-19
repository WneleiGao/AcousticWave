function OneStepAdjoint!(spt2::SnapShot, spt1::SnapShot, fidMtx::FidMtx)
    spt2.it = spt1.it - 1
    spt2.vz[:] = spt1.vz + (fidMtx.MpzBvz)' * spt1.pz
    spt2.vx[:] = spt1.vx + (fidMtx.MpxBvx)' * spt1.px
    spt2.pz[:] =           (fidMtx.MpzBpz)' * spt1.pz
    spt2.px[:] =           (fidMtx.MpxBpx)' * spt1.px

    spt2.pz[:] = (fidMtx.MvzBp )' * spt2.vz + (fidMtx.MvxBp)' * spt2.vx + spt2.pz
    spt2.px[:] = (fidMtx.MvzBp )' * spt2.vz + (fidMtx.MvxBp)' * spt2.vx + spt2.px
    spt2.vz[:] = (fidMtx.MvzBvz)' * spt2.vz
    spt2.vx[:] = (fidMtx.MvxBvx)' * spt2.vx
    return nothing
end

# save part(ntw ~ 1) of adjoint spt, vp, p, input shot
function MultiStepAdjoint(path::String , ntw::Int64, shot::Shot, fidMtx::FidMtx, wtype="p")
    nz = fidMtx.nz; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot  ; dt = shot.dt  ; nt = shot.nt   ;
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot)
    path_tmp = join([path "_tmp"])
    if ntw > nt
       error("length of src is larger than records")
    elseif ntw == nt
       if wtype == "spt"
          fid = writeSnapShot(path_tmp, spt1)
       elseif wtype == "vp"
          fid = writeWfd(path_tmp, spt1)
       elseif wtype == "p"
          fid = writeStress(path_tmp, spt1)
       end
       for it = ntw-1 : -1 : 1
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot2SnapShot!(spt2, shot)
           copySnapShot!(spt1, spt2)
           if wtype == "spt"
              writeSnapShot(fid, spt1)
           elseif wtype == "vp"
              writeWfd(fid, spt1)
           elseif wtype == "p"
              writeStress(fid, spt1)
           end
       end
    elseif ntw < nt
       for it = nt-1: -1: ntw
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot2SnapShot!(spt2, shot)
           copySnapShot!(spt1, spt2)
       end
       if wtype == "spt"
          fid = writeSnapShot(path_tmp, spt1)
       elseif wtype == "vp"
          fid = writeWfd(path_tmp, spt1)
       elseif wtype == "p"
          fid = writeStress(path_tmp, spt1)
       end
       for it = ntw-1 : -1 : 1
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot2SnapShot!(spt2, shot)
           copySnapShot!(spt1, spt2)
           if wtype == "spt"
              writeSnapShot(fid, spt1)
           elseif wtype == "vp"
              writeWfd(fid, spt1)
           elseif wtype == "p"
              writeStress(fid, spt1)
           end
       end
    end
    close(fid)
    if wtype == "spt"
       reverseSnapShotsOrder(path, path_tmp)
    elseif wtype == "vp"
       reverseWfdOrder(path, path_tmp)
    elseif wtype == "p"
       reverseStressOrder(path, path_tmp)
    end
    return nothing
end

# save part(ntw ~ 1) of adjoint spt, vp, p, input shot3
function MultiStepAdjoint(path::String , ntw::Int64, shot::Shot3, fidMtx::FidMtx; wtype="p")
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; dt = shot.dt  ; nt = shot.nt;
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot32SnapShot!(spt1, shot)
    path_tmp = join([path "_tmp"])
    if ntw > nt
       error("length of src is larger than records")
    elseif ntw == nt
       if wtype == "spt"
          fid = writeSnapShot(path_tmp, spt1)
       elseif wtype == "vp"
          fid = writeWfd(path_tmp, spt1)
       elseif wtype == "p"
          fid = writeStress(path_tmp, spt1)
       end
       for it = ntw-1 : -1 : 1
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot32SnapShot!(spt2, shot)
           copySnapShot!(spt1, spt2)
           if wtype == "spt"
              writeSnapShot(fid, spt1)
           elseif wtype == "vp"
              writeWfd(fid, spt1)
           elseif wtype == "p"
              writeStress(fid, spt1)
           end
       end
    elseif ntw < nt
       for it = nt-1: -1: ntw
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot32SnapShot!(spt2, shot)
           copySnapShot!(spt1, spt2)
       end
       if wtype == "spt"
          fid = writeSnapShot(path_tmp, spt1)
       elseif wtype == "vp"
          fid = writeWfd(path_tmp, spt1)
       elseif wtype == "p"
          fid = writeStress(path_tmp, spt1)
       end
       for it = ntw-1 : -1 : 1
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot32SnapShot!(spt2, shot)
           copySnapShot!(spt1, spt2)
           if wtype == "spt"
              writeSnapShot(fid, spt1)
           elseif wtype == "vp"
              writeWfd(fid, spt1)
           elseif wtype == "p"
              writeStress(fid, spt1)
           end
       end
    end
    close(fid)
    if wtype == "spt"
       reverseSnapShotsOrder(path, path_tmp)
    elseif wtype == "vp"
       reverseWfdOrder(path, path_tmp)
    elseif wtype == "p"
       reverseStressOrder(path, path_tmp)
    end
    return nothing
end

# output image I, input shot and source side wavefield (only stress)
function MultiStepAdjoint(shot::Shot, path::String , fidMtx::FidMtx)
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; dt = shot.dt  ; nt = shot.nt;
    (nz1, nx1, ext1, iflag1, dt1, nt1) = InfoStress(path)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot)
    I = zeros(nz, nx)
    if nt == nt1
       I = RTMimage(I, spt1, path)
    end
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot2SnapShot!(spt2, shot)
        copySnapShot!(spt1, spt2)
        if it <= nt1
           I = RTMimage(I, spt1, path)
        end
    end
    return I
end

# output image I, input shot3 and source side wavefield (only stress)
function MultiStepAdjoint(shot::Shot3, path::String , fidMtx::FidMtx)
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; dt = shot.dt  ; nt = shot.nt;
    (nz1, nx1, ext1, iflag1, dt1, nt1) = InfoStress(path)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot32SnapShot!(spt1, shot)
    I = zeros(nz, nx)
    if nt == nt1
       I = RTMimage(I, spt1, path)
    end
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot32SnapShot!(spt2, shot)
        copySnapShot!(spt1, spt2)
        if it <= nt1
           I = RTMimage(I, spt1, path)
        end
    end
    return I
end

# output the distribution of sources, assume the wavelet is known, input shot
function MultiStepAdjoint(w::Array{Float64,1}, shot::Shot, fidMtx::FidMtx)
    nz = fidMtx.nz; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot  ; dt = shot.dt  ; nt = shot.nt   ; ntw   = length(w)   ;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    dis= zeros(nz, nx)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot)
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot2SnapShot!(spt2, shot)
        copySnapShot!(spt1, spt2)
        if it <= ntw
           dis = dis + spt2dis(w, spt1)
        end
    end
    return dis
end

# output the distribution of sources, input shot3, assume the wavelet is known
function MultiStepAdjoint(w::Array{Float64,1}, shot::Shot3, fidMtx::FidMtx)
    nz = fidMtx.nz; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot  ; dt = shot.dt  ; nt = shot.nt   ; ntw   = length(w)   ;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    dis= zeros(nz, nx)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot32SnapShot!(spt1, shot)
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot32SnapShot!(spt2, shot)
        copySnapShot!(spt1, spt2)
        if it <= ntw
           dis = dis + spt2dis(w, spt1)
        end
    end
    return dis
end

# out the wavelet, assume the distribution of source is known, input shot
function MultiStepAdjoint(ntw::Int64, dis::Array{Float64,2}, shot::Shot, fidMtx::FidMtx)
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; dt = shot.dt  ; nt = shot.nt   ;
    if iflag == 1
       zupper = ext
       Nz = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    w  = zeros(ntw)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot)
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot2SnapShot!(spt2, shot)
        copySnapShot!(spt1, spt2)
        if it <= ntw
           spt2wlet(w, dis, spt1)
        end
    end
    return w
end

# out the wavelet, assume the distribution of source is known, input shot3
function MultiStepAdjoint(ntw::Int64, dis::Array{Float64,2}, shot::Shot3, fidMtx::FidMtx)
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; dt = shot.dt  ; nt = shot.nt   ;
    if iflag == 1
       zupper = ext
       Nz = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    w  = zeros(ntw)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot32SnapShot!(spt1, shot)
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot32SnapShot!(spt2, shot)
        copySnapShot!(spt1, spt2)
        if it <= ntw
           spt2wlet(w, dis, spt1)
        end
    end
    return w
end

function MultiStepAdjoint(oz::Int64, ox::Int64, wnz::Int64, wnx::Int64, wnt::Int64, shot::Shot, fidMtx::FidMtx)
    nz = fidMtx.nz; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag;
    dt = shot.dt;   nt = shot.nt  ;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    wsc = InitWSrcCube(oz, ox, wnz, wnx, ext, iflag, 0.0, dt, wnt)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot)
    for it = nt-1 : -1 : 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot2SnapShot!(spt2, shot)
        copySnapShot!(spt1, spt2)
        if it <= wnt
           spt2Wsc!(wsc, spt1)
        end
    end
    return wsc
end

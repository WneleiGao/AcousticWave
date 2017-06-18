using PyPlot, AcousticWave
nz = 300; nx = 300; ext = 20; iflag = 1;
dz = 15.; dx = 15.; dt = 1.e-3; f0=10.0; tmax = 10.0;
v  = 2500. * ones(nz, nx);
v[250:end,:]  = 3500.
# initialize source term
isz = 5; isx = 250; ot = 0.; amp = 1.0;
src = InitSource(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# discretize spatial derivative operator
vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, v, f0);
irx = collect(1:2:nx); irz = 5*ones(Int64, length(irx));
@time shot = MultiStepForward(irz, irx, src, fidMtx, tmax=tmax);
@time shot1 = MultiStepForward_test(irz, irx, src, fidMtx, tmax=tmax)


function MultiStepForward_test(irz::Array{Int64,1}, irx::Array{Int64,1}, srcs::Array{Source,1}, fidMtx::FidMtx; tmax=1.0)
    nz = fidMtx.nz ; nx = fidMtx.nx; ext= fidMtx.ext; iflag = fidMtx.iflag; dt = fidMtx.dt;
    nt = round(Int64, tmax/dt)+1
    (tl, tu) = SrcRange(srcs)
    shot = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt)
    spt1 = InitSnapShot(nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(nz, nx, ext, iflag, dt, 2)
    AddMultiSources!(spt1, srcs)
    spt2shot!(shot, spt1)
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx)
        if tl <= (it-1)*dt <= tu
           AddMultiSources!(spt2, srcs)
        end
        copySnapShot!(spt1, spt2)
        spt2shot!(shot, spt1)
    end
    return shot
end


# spt1 = InitSnapShot(nz, nx, ext, 1, dt, 1);
# spt1.vz = randn(length(spt1.vz));
# spt1.vx = randn(length(spt1.vx));
# spt1.pz = randn(length(spt1.pz));
# spt1.px = randn(length(spt1.px));
# spt2 = InitSnapShot(nz, nx, ext, 1, dt, 1);
# @time for i = 1 : 1429
#     OneStepForward!(spt2, spt1, fidMtx);
# end
# tmp = rand(length(spt1.vz));
# tmp1 = rand(length(spt1.vz));
# @time for i = 1 : 1429
#      OneStepForward!(spt2, spt1, fidMtx, tmp, tmp1);
# end







root = homedir(); path = join([root "/Desktop/stress.bin"])
MultiStepForward(path, src, fidMtx, tmax=tmax)
shot = stress2Shot(irz, irx, path)
SeisPlot(shot.p, pclip=95, cmap="PuOr", oy=dt, ox=dx)

path1 = join([root "/Desktop/stress"])
waveAnim(path, path1, dx=dx, dz=dz, interval=50)


# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 20  ; iflag = 1 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx);
# v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125, :] = 4000.;
# v[126:nz, :]  = 4500.;
#
# vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# # finish testing forward modeling code
# isx=collect(30:30:180); isz=ones(Int64,length(isx)); ot=zeros(length(isx)); amp=ones(length(isx));
# srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# root = homedir(); path = join([root "/Desktop/srcp.bin"])
# writeSrcs2Stress(path, srcs)
# irx = collect(1:3:nx); irz = ones(Int64, length(irx));
# shot3 = MultiStepForward(irz, irx, path, fidMtx, tmax=tmax, otype="vp");
# shotc = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax, otype="vp");
# dis = zeros(nz,nx); dis[isz,isx] = 1.0; w = Ricker(f0, dt);
# shotb = MultiStepForward(irz, irx, w, dis, fidMtx, tmax=tmax, otype="vp")
# vecnorm(shot3.vz-shotc.vz)
# vecnorm(shot3.vx-shotc.vx)
# vecnorm(shot3.p -shotc.p )
# vecnorm(shotb.vz-shotc.vz)
# vecnorm(shotb.vx-shotc.vx)
# vecnorm(shotb.p -shotc.p )
# # I = ones(nz, nx)*1/2;
# # shotb = MultiStepForward_test(irz, irx, I, path, fidMtx, tmax=tmax, otype="p");
# vecnorm(shotb.vz*2-shotc.vz)
# vecnorm(shotb.vx*2-shotc.vx)
# vecnorm(shotb.p*2 -shotc.p )
# val = quantile(abs(vec(shotb.p)), 0.90);
# imshow(shotb.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.04);
# =================================================================
# output shot, input single source
# isz = 1; isx = 100; ot = 0.0; amp = 1.0;
# src = InitSource(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# irx = collect(1:2:nx); irz = ones(Int64, length(irx));
# shot3 = MultiStepForward(irz, irx, src, fidMtx, tmax=tmax, otype="vp");
# val = quantile(abs(vec(shot3.p)), 0.95);
# imshow(shot3.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.04);
# root = homedir(); path = join([root "/Desktop/acoustic_wfd.bin"])
# MultiStepForward(path, src, fidMtx, tmax=tmax, wtype="wfd");
# shotc = wfd2Shot3_test(irz, irx, path);
# vecnorm(shotc.vz-shot3.vz);
# vecnorm(shotc.vx-shot3.vx);
# vecnorm(shotc.p-shot3.p);
# val = quantile(abs(vec(shot3.vz)), 0.95);
# imshow(shot3.vz, vmin=-val, vmax=val, cmap="PuOr", aspect=0.04);
# val = quantile(abs(vec(shot3.vx)), 0.95);
# imshow(shot3.vx, vmin=-val, vmax=val, cmap="PuOr", aspect=0.04);
# val = quantile(abs(vec(shot3.p)), 0.95);
# imshow(shotc.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.04);

# =================================================================
# save data to slow memory, input multiple sources
# isx=collect(50:50:150); isz=ones(Int64,length(isx)); ot=zeros(length(isx)); amp=ones(length(isx));
# srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# root = homedir(); path = join([root "/Desktop/acoustic.bin"])
# MultiStepForward(path, srcs, fidMtx, tmax=tmax, wtype="spt");
# root = homedir(); path = join([root "/Desktop/acoustic_vp.bin"])
# MultiStepForward(path, srcs, fidMtx, tmax=tmax, wtype="vp");
# root = homedir(); path = join([root "/Desktop/acoustic_p.bin"])
# MultiStepForward(path, srcs, fidMtx, tmax=tmax, wtype="p");
# it = 347;
# path=join([root "/Desktop/acoustic.bin"]); spt=readSnapShot(path,it);
# path=join([root "/Desktop/acoustic_vp.bin"]); wfd=readWfd(path,it);
# path=join([root "/Desktop/acoustic_p.bin"]);  sts=readStress(path,it);
# Nz = nz+ext; Nx = nx + 2*ext
# vz = reshape(spt.vz, Nz, Nx)[1:nz, ext+1:ext+nx];
# vx = reshape(spt.vx, Nz, Nx)[1:nz, ext+1:ext+nx];
# p = reshape(spt.pz+spt.px, Nz, Nx)[1:nz, ext+1:ext+nx];
# vecnorm(p-sts.p)
# vecnorm(vz-wfd.vz)
# vecnorm(vx-wfd.vx)
# vecnorm(wfd.p-sts.p)

# ==============================================================================
# dot product test for one step
# spt1 = InitSnapShot(101, 101, nz, nx, ext, 1, dt, 1);
# spt1.vz = randn(length(spt1.vz));
# spt1.vx = randn(length(spt1.vx));
# spt1.pz = randn(length(spt1.pz));
# spt1.px = randn(length(spt1.px));
# spt2 = InitSnapShot(101, 101, nz, nx, ext, 1, dt, 1);
# OneStepForward!(spt2, spt1, fidMtx);
#
# tmp2 = InitSnapShot(101, 101, nz, nx, ext, 1, dt, 2)
# tmp2.vz = randn(length(tmp2.vz))
# tmp2.vx = randn(length(tmp2.vx))
# tmp2.pz = randn(length(tmp2.pz))
# tmp2.px = randn(length(tmp2.px))
# tmp1 = InitSnapShot(101, 101, nz, nx, ext, 1, dt, 2)
# OneStepAdjoint!(tmp1, tmp2, fidMtx)
#
# x  = vcat(spt1.vz, spt1.vx, spt1.pz, spt1.px)
# x1 = vcat(tmp1.vz, tmp1.vx, tmp1.pz, tmp1.px)
# y  = vcat(spt2.vz, spt2.vx, spt2.pz, spt2.px)
# y1 = vcat(tmp2.vz, tmp2.vx, tmp2.pz, tmp2.px)
#
# x'*x1 - y'*y1

# ==============================================================================
# nz = 201; nx = 201; dt = 2e-3; iflag = 1; ext = 20
# if iflag ==1
#    zupper = ext
#    Nz = nz + 2*ext
# elseif iflag == 2
#    Nz = nx +   ext
#    zupper = ext
# end
# Nx = nx + 2*ext
# spt1 = InitSnapShot(101, 101, nz, nx, ext, iflag, dt, 1);
# spt1.vz = randn(length(spt1.vz));
# spt1.vx = randn(length(spt1.vx));
# spt1.pz = randn(length(spt1.pz));
# spt1.px = randn(length(spt1.px));
#
# iz = collect(1:5:201); ix = collect(1:3:201);
# nrz= length(iz); nrx = length(ix);
# iz = vec(repmat(iz, nrx, 1));
# ix = vec(repmat(ix', (nrz), 1));
# pos = hcat(iz, ix);
# shot1 = InitShot(101, 101, pos, 0.0, 1, dt)
# spt2shot!(shot1, spt1, Nz, Nx, ext, zupper)
#
#
# shot2 = InitShot(101, 101, pos, 0.0, 1, dt)
# shot2.d = randn(1, size(shot2.d, 2))
# spt2 = InitSnapShot(101, 101, nz, nx, ext, iflag, dt, 1);
# AddShot2SnapShot!(spt2, shot2, Nz, Nx, ext, zupper)
#
# x  = vcat(spt1.vz, spt1.vx, spt1.pz, spt1.px)
# x1 = vcat(spt2.vz, spt2.vx, spt2.pz, spt2.px)
# y  = vec(shot1.d)
# y1 = vec(shot2.d)
# x'*x1 - y'*y1


# ==============================================================================
# using PyPlot, AcousticWave
# nz = 200; nx = 200; ext = 20; iflag = 2;
# dz = 15.; dx = 15.; dt = 1.5e-3; f0=10.0; tmax = 4.0;
# v  = 2000. * ones(nz, nx);
# v[101:end, :] = 3000.
# # initialize source term
# isz = 1; isx = 120; ot = 0.0; amp = 1.0;
# src = InitSource(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# # discretize spatial derivative operator
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, v, f0);
# irx = collect(1:2:nx); irz = ones(Int64, length(irx));
# shot3 = MultiStepForward(irz, irx, src, fidMtx, tmax=tmax, otype="vp");
# SeisPlot(shot3.p, pclip=90)
# path = "/Users/wenlei/Desktop/acoustic.bin"
# MultiStepForward(path, src, fidMtx, tmax=1.0);

# iz = 5 * ones(Int64, 201); ix = collect(1:201);
# pos = hcat(iz, ix)
# shot = fullWfd2Shot(pos, path);
# shot1 = MultiStepForward(pos, src[1], fidMtx, tmax=1.0)
#
# path = "/Users/wenlei/Desktop/adjoint.bin"
# MultiStepAdjoint!(path, shot, fidMtx)

# ==============================================================================
# using AcousticWave
# nz = 201; nx = 201; ext = 20; iflag = 1; dz = 10.; dx = 10.; dt = 2e-3; f0 = 30.0
# v  = 3000. * ones(nz, nx);
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# # initialize source term
# ix  = collect(1:2:201); iz = 101*ones(Int64, length(ix));
# pos = hcat(iz, ix)
# f0 = 30.0; ot = zeros(size(pos,1));
# src = InitSources(pos, f0, ot, dt);
# path = "/Users/wenlei/Desktop/src2spts.bin"
# src2SnapShots(path, src, nz, nx, ext, iflag);
#
# # position of receiver
# irz = []; irx = []
# for ix = 1 : nx
#     for iz = 1 : nz
#         if rand() <= 0.1
#            irz = push!(irz, iz)
#            irx = push!(irx, ix)
#         end
#     end
# end
# pos_rec = convert(Array{Int64,2}, hcat(irz, irx))
# shot  = MultiStepForward(pos_rec, path, fidMtx, tmax=1.0)
# path_wfd = "/Users/wenlei/Desktop/acoustic.bin"
# MultiStepForward(path_wfd, src , fidMtx, tmax=1.0)
# shot1 = fullWfd2Shot(pos_rec, path_wfd)

# ==============================================================================
# using AcousticWave
# nz = 201; nx = 201; ext = 20; iflag = 1; dz = 10.; dx = 10.; dt = 2e-3; f0 = 30.0
# v  = 3000. * ones(nz, nx);
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# # initialize source term
# isz = []; isx = []
# for ix = 1 : nx
#     for iz = 1 : nz
#         if rand() <= 0.1
#            isz = push!(isz, iz)
#            isx = push!(isx, ix)
#         end
#     end
# end
# pos_src = convert(Array{Int64,2}, hcat(isz, isx))
# f0 = 30.0; ot = zeros(size(pos_src,1));
# src = InitSources(pos_src, f0, ot, dt);
# tmax = 1; nt = round(Int64, tmax/dt+1);
# w = Ricker(f0, dt); ns = nt + 1 - length(w)
# for isrc = 1 : length(src)
#     src[isrc].nt= nt
#     src[isrc].p = conv(randn(ns),w)
# end
# path = "/Users/wenlei/Desktop/src2spts.bin"
# src2SnapShots(path, src, nz, nx, ext, iflag);
# # position of receiver
# irx = collect(1:2:nx); irz = 5*ones(Int64,length(irx));
# pos_rec = hcat(irz, irx)
# shot  = MultiStepForward(pos_rec, path, fidMtx, tmax=1.0)
#
# path_adj = "/Users/wenlei/Desktop/adjoint.bin"
# shot1 = InitShot(0, 0, pos_rec, 0.0, nt, dt)
# for irec = 1 : length(shot1.irz)
#     shot1.d[:,irec] = conv(randn(ns), w)
# end
# MultiStepAdjoint_test(path_adj, shot1, fidMtx)
#
# L = 0.0
# for ir = 1 : length(shot.irx)
#     L = L + dot(shot.d[:,ir], shot1.d[:,ir])
# end
#
# T = 0.0
# for it = 1 : nt
#     spt1 = readSnapShot(path    , it)
#     spt2 = readSnapShot(path_adj, it)
#     T = T + dot(spt1.vz, spt2.vz) + dot(spt1.vx, spt2.vx) + dot(spt1.pz, spt2.pz) + dot(spt1.px, spt2.px)
# end

# ==============================================================================
# using AcousticWave
# path = "/Users/wenlei/Desktop/acoustic.bin"
# nz = 201; nx = 201; ext = 20; iflag = 1; dz = 10.; dx = 10.; dt = 2e-3;
# v  = 3000. * ones(nz, nx);
# # initialize source term
# pos = [101 101]; f0 = 30.0; ot = [0.0];
# src = InitSources(pos, f0, ot, dt);
# # discretize spatial derivative operator
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# path = "/Users/wenlei/Desktop/acoustic.bin"
# MultiStepForward(path, src, fidMtx, tmax=1.0);
#
# using PyPlot, AcousticWave
# path = "/Users/wenlei/Desktop/acoustic.bin"
# pathout = "/Users/wenlei/Desktop/homo"
# waveAnim(path, pathout)
# sptPlot(path, 100)
# irx = collect(1:2:nx); irz = 5*ones(Int64,length(irx));
# pos_rec = hcat(irz, irx)
# shot = fullWfd2Shot(pos_rec, path)

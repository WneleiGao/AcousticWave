# ======dot product test, input stress, adjoint is also stress=======
using PyPlot, AcousticWave
nz = 150; nx = 200; ext= 30  ; iflag = 1 ;
dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
for iz = 1 : 25
    v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
end
v[101:125, :] = 4000.; v[126:nz, :] = 4500.;

vl = minimum(v); vu = maximum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# initialize source term
isz=[]; isx=[]; ot=[]; amp=[];
for ix = 1 : nx
    for iz = 1 : nz
        if rand() <= 0.05
           push!(isz, iz)
           push!(isx, ix)
           push!(ot, floor(Int64,rand()/dt)*dt)
           if rand() >= 0.5
              push!(amp, 1.0)
           else
              push!(amp,-1.0)
           end
        end
    end
end
isz = convert(Array{Int64},isz); isx = convert(Array{Int64},isx);
ot  = convert(Array{Float64},ot ); amp = convert(Array{Float64},amp);
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
(tl, tu) = SrcRange(srcs); ntw = round(Int64, tu/dt)+1;
root = homedir(); path = join([root "/Desktop/srcp.bin"]);
writeSrcs2Stress(path, srcs)

irx = collect(1:nx); irz = ones(Int64, length(irx));
shot = MultiStepForward(irz, irx, path, fidMtx, tmax=tmax);
val = quantile(abs(vec(shot.p)), 0.90); imshow(shot.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.1);

root = homedir(); path1 = join([root "/Desktop/adj.bin"]);
nt = size(shot.p ,1);
shot1 = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
wlet=Ricker(f0,dt); nw=length(wlet);
for ir = 1 : length(irz)
    if rand() < 0.7
       if rand() >= 0.5
          tmp = 1.0 * wlet
       else
          tmp =-1.0 * wlet
       end
       inds = floor(Int64, rand()*(nt-nw))
       shot1.p[inds:inds+nw-1 ,ir] = tmp[:]
    end
end
imshow(shot1.p, aspect=0.1);
MultiStepAdjoint(path1, ntw, shot1, fidMtx)

L = 0.0
for ir = 1 : length(shot.irx)
    L = L + dot(shot.p[:,ir], shot1.p[:,ir])
end

T = 0.0
for it = 1 : ntw
    spt1 = readStress(path , it)
    spt2 = readStress(path1, it)
    T = T + dot(vec(spt1.p), vec(spt2.p)*1/2) # determined by the source adding style
end

# ======dot product test, forward produce Shot3, adjoint produce stress=======
# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 30  ; iflag = 2 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125, :] = 4000.; v[126:nz, :] = 4500.;
# vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
#
# # initialize source term
# isz=[]; isx=[]; ot=[]; amp=[];
# for ix = 1 : nx
#     for iz = 1 : nz
#         if rand() <= 0.05
#            push!(isz, iz)
#            push!(isx, ix)
#            push!(ot, floor(Int64,rand()/dt)*dt)
#            if rand() >= 0.5
#               push!(amp, 1.0)
#            else
#               push!(amp,-1.0)
#            end
#         end
#     end
# end
# isz = convert(Array{Int64},isz);   isx = convert(Array{Int64},isx)  ;
# ot  = convert(Array{Float64},ot ); amp = convert(Array{Float64},amp);
# srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# (tl, tu) = SrcRange(srcs); ntw = round(Int64, tu/dt)+1;
# root = homedir(); path = join([root "/Desktop/srcp.bin"]);
# writeSrcs2Stress(path, srcs);
#
# irx  = collect(1:nx); irz = ones(Int64, length(irx));
# shot3 = MultiStepForward(irz, irx, path, fidMtx, tmax=tmax, otype="vp");
# val = quantile(abs(vec(shot3.vx)), 0.90); imshow(shot3.vx, vmin=-val, vmax=val, cmap="PuOr", aspect=0.1);
#
# root = homedir(); path1 = join([root "/Desktop/adj.bin"]);
# nt = size(shot3.p ,1);
# shot3_adj = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
# wlet=Ricker(f0,dt); nw=length(wlet);
# for ir = 1 : length(irz)
#     if rand() < 0.7
#        if rand() >= 0.5
#           tmp = 1.0 * wlet
#        else
#           tmp =-1.0 * wlet
#        end
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.vz[inds:inds+nw-1 ,ir] = tmp[:]
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.vx[inds:inds+nw-1 ,ir] = tmp[:]
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.p[inds:inds+nw-1 ,ir] = tmp[:]
#     end
# end
# imshow(shot3_adj.p, aspect=0.1); figure(); imshow(shot3_adj.vz, aspect=0.1); figure(); imshow(shot3_adj.vx, aspect=0.1);
# MultiStepAdjoint(path1, ntw, shot3_adj, fidMtx);
#
# L = 0.0
# for ir = 1 : length(shot3.irx)
#     L = L + dot(shot3.p[:,ir],shot3_adj.p[:,ir]) + dot(shot3.vz[:,ir],shot3_adj.vz[:,ir]) + dot(shot3.vx[:,ir],shot3_adj.vx[:,ir])
# end
#
# T = 0.0
# for it = 1 : ntw
#     spt1 = readStress(path , it)
#     spt2 = readStress(path1, it)
#     T = T + dot(vec(spt1.p), vec(spt2.p)*1/2) # determined by the source adding style
# end

# ======dot product test, forward produce Shot3, adjoint produce stress=======
# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 30  ; iflag = 2 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125,:] = 4000.; v[126:nz,:] = 4500.; vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
#
# isx=collect(1:nx); isz=ones(Int64,length(isx)); ot=zeros(length(isx)); amp=ones(length(isx));
# srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# root = homedir(); path = join([root "/Desktop/srcwfd.bin"])
# MultiStepForward(path, srcs, fidMtx, tmax=tmax)
#
# irx = collect(1:nx); irz = ones(Int64, length(irx)); I = randn(nz, nx);
# shot3 = MultiStepForward(irz, irx, I, path, fidMtx, tmax=tmax, otype="vp");
# val = quantile(abs(vec(shot3.vx)), 0.90); imshow(shot3.vx, vmin=-val, vmax=val, cmap="PuOr", aspect=0.1);
#
# # ===adjoint===
# nt = size(shot3.p ,1);
# shot3_adj = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
# wlet=Ricker(f0,dt); nw=length(wlet);
# for ir = 1 : length(irz)
#     if rand() < 0.7
#        if rand() >= 0.5
#           tmp = 1.0 * wlet
#        else
#           tmp =-1.0 * wlet
#        end
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.vz[inds:inds+nw-1 ,ir] = tmp[:]
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.vx[inds:inds+nw-1 ,ir] = tmp[:]
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.p[inds:inds+nw-1 ,ir] = tmp[:]
#     end
# end
# Iadj = MultiStepAdjoint(shot3_adj, path, fidMtx);
# L = 0.0
# for ir = 1 : length(shot3.irx)
#     L = L + dot(shot3.p[:,ir],shot3_adj.p[:,ir]) + dot(shot3.vz[:,ir],shot3_adj.vz[:,ir]) + dot(shot3.vx[:,ir],shot3_adj.vx[:,ir])
# end
# T = dot(vec(I), vec(Iadj));

#=====================dot_product test for imaging, recording is shot==============================#
# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 30  ; iflag = 1 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125, :] = 4000.; v[126:nz, :] = 4500.;
# vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# # produce source side wave field
# isx=collect(1:nx); isz=ones(Int64,length(isx)); ot=zeros(length(isx)); amp=ones(length(isx));
# srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# root = homedir(); path = join([root "/Desktop/srcwfd.bin"])
# MultiStepForward(path, srcs, fidMtx, tmax=tmax)
#
# I = randn(nz, nx); irx = collect(1:nx); irz = ones(Int64, length(irx));
# shot = MultiStepForward(irz, irx, I, path, fidMtx, tmax=tmax);
# # val = quantile(abs(vec(shot.p)), 0.90); imshow(shot.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.1);
#
# nt = size(shot.p ,1);
# shot1 = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
# wlet=Ricker(f0,dt); nw=length(wlet);
# for ir = 1 : length(irz)
#     if rand() < 0.7
#        if rand() >= 0.5
#           tmp = 1.0 * wlet
#        else
#           tmp =-1.0 * wlet
#        end
#        inds = floor(Int64, rand()*(nt-nw))
#        shot1.p[inds:inds+nw-1 ,ir] = tmp[:]
#     end
# end
# I1 = MultiStepAdjoint(shot1, path, fidMtx)
# L = dot(vec(I), vec(I1))
# T = dot(vec(shot.p), vec(shot1.p))


# ===========adjoint_test, input shot3 and shot, produce same result===================================
# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 30  ; iflag = 2 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125, :] = 4000.; v[126:nz, :] = 4500.;
# vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# # initialize source term
# isz=[]; isx=[]; ot=[]; amp=[];
# for ix = 1 : nx
#     for iz = 1 : nz
#         if rand() <= 0.05
#            push!(isz, iz)
#            push!(isx, ix)
#            push!(ot, floor(Int64,rand()/dt)*dt)
#            if rand() >= 0.5
#               push!(amp, 1.0)
#            else
#               push!(amp,-1.0)
#            end
#         end
#     end
# end
# isz = convert(Array{Int64},isz); isx = convert(Array{Int64},isx);
# ot  = convert(Array{Float64},ot ); amp = convert(Array{Float64},amp);
# srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
# (tl, tu) = SrcRange(srcs); ntw = round(Int64, tu/dt)+1;
# irx  = collect(1:nx); irz = ones(Int64, length(irx));
# shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax);
# val = quantile(abs(vec(shot.p)), 0.90); imshow(shot.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.1);
# root = homedir(); path1 = join([root "/Desktop/adj1.bin"]);
# MultiStepAdjoint(path1, ntw, shot, fidMtx);
#
# root = homedir(); path2 = join([root "/Desktop/adj2.bin"]);
# nt = size(shot.p ,1);
# shot3 = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
# shot3.p[:] = shot.p[:]
# MultiStepAdjoint(path2, ntw, shot3, fidMtx);
# T = 0.0
# for it = 1 : ntw
#     spt1 = readStress(path1, it)
#     spt2 = readStress(path2, it)
#     tmp  = vec(spt1.p) - vec(spt2.p)
#     T = T + dot(tmp, tmp) # determined by the source adding style
# end

# # ============dpt, forward produce shot, input distribution of source=======================
# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 30  ; iflag = 1 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125, :] = 4000.; v[126:nz, :] = 4500.;
# vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# # produce source side wave field
# dis = randn(nz,nx); w = Ricker(f0, dt);
# irx = collect(1:nx); irz = ones(Int64, length(irx));
# shot = MultiStepForward(irz, irx, w, dis, fidMtx, tmax=tmax);
# # val = quantile(abs(vec(shot.p)), 0.90); imshow(shot.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.1);
#
# nt = size(shot.p ,1); nw = length(w)
# shot_adj = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
# for ir = 1 : length(irz)
#     if rand() < 0.7
#        if rand() >= 0.5
#           tmp = 1.0 * w
#        else
#           tmp =-1.0 * w
#        end
#        inds = floor(Int64, rand()*(nt-nw))
#        shot_adj.p[inds:inds+nw-1 ,ir] = tmp[:]
#     end
# end
# dis_adj = MultiStepAdjoint(w, shot_adj, fidMtx)
# L = dot(vec(shot.p), vec(shot_adj.p))
# T = dot(vec(dis), vec(dis_adj))

# =====dpt, forward produce shot3, input distribution of source (not perfect)=====
# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 30  ; iflag = 1 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125, :] = 4000.; v[126:nz, :] = 4500.;
# vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# # produce source side wave field
# dis = randn(nz,nx); w = Ricker(f0, dt);
# irx = collect(1:nx); irz = ones(Int64, length(irx));
# shot3 = MultiStepForward(irz, irx, w, dis, fidMtx, tmax=tmax, otype="vp");
# # val = quantile(abs(vec(shot.p)), 0.90); imshow(shot.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.1);
# nt = size(shot3.p ,1); nw = length(w);
# shot3_adj = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
# for ir = 1 : length(irz)
#     if rand() < 0.7
#        if rand() >= 0.5
#           tmp = 1.0 * w
#        else
#           tmp =-1.0 * w
#        end
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.vz[inds:inds+nw-1 ,ir] = tmp[:]
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.vx[inds:inds+nw-1 ,ir] = tmp[:]
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot3_adj.p[inds:inds+nw-1 ,ir] = tmp[:]
#     end
# end
# dis_adj = MultiStepAdjoint(w, shot3_adj, fidMtx)
# L = dot(vec(shot3.p), vec(shot3_adj.p)) + dot(vec(shot3.vz), vec(shot3_adj.vz)) + dot(vec(shot3.vx), vec(shot3_adj.vx))
# T = dot(vec(dis), vec(dis_adj))


# ====dpt, forward produce shot, input source wavelet, assume the distribution is known====
# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 30  ; iflag = 1 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125, :] = 4000.; v[126:nz, :] = 4500.;
# vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# # produce source side wave field
# dis = randn(nz,nx); w = Ricker(f0, dt);
# irx = collect(1:nx); irz = ones(Int64, length(irx));
# shot = MultiStepForward(irz, irx, w, dis, fidMtx, tmax=tmax);
#
# nt = size(shot.p, 1); nw = length(w);
# shot_adj = InitShot(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
# for ir = 1 : length(irz)
#     if rand() < 0.7
#        if rand() >= 0.5
#           tmp = 1.0 * w
#        else
#           tmp =-1.0 * w
#        end
#        inds = floor(Int64, rand()*(nt-nw-100))
#        shot_adj.p[inds:inds+nw-1 ,ir] = tmp[:]
#     end
# end
# wadj = MultiStepAdjoint(nw, dis, shot_adj, fidMtx)
# L = dot(vec(shot.p), vec(shot_adj.p))
# T = dot(w, wadj)

# ====dpt, forward produce shot, input source wavelet, assume the distribution is known====
# using PyPlot, AcousticWave
# nz = 150; nx = 200; ext= 30  ; iflag = 1 ;
# dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
# v = 2500.*ones(nz, nx); v[26:100, :] = 3000.;
# for iz = 1 : 25
#     v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
# end
# v[101:125, :] = 4000.; v[126:nz, :] = 4500.;
# vl = minimum(v); vu = maximum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# # produce source side wave field
# dis = randn(nz,nx); w = Ricker(f0, dt);
# irx = collect(1:nx); irz = ones(Int64, length(irx));
# shot = MultiStepForward(irz, irx, w, dis, fidMtx, tmax=tmax, otype="vp");
#
# nt = size(shot.p, 1); nw = length(w);
# shot_adj = InitShot3(0, 0, nz, nx, ext, iflag, irz, irx, 0.0, dt, nt);
# for ir = 1 : length(irz)
#     if rand() < 0.7
#        if rand() >= 0.5
#           tmp = 1.0 * w
#        else
#           tmp =-1.0 * w
#        end
#        inds = floor(Int64, rand()*(nt-nw-150))
#        shot_adj.vz[inds:inds+nw-1 ,ir] = tmp[:]
#        inds = floor(Int64, rand()*(nt-nw-150))
#        shot_adj.vx[inds:inds+nw-1 ,ir] = tmp[:]
#        inds = floor(Int64, rand()*(nt-nw-150))
#        shot_adj.p[inds:inds+nw-1 ,ir] = tmp[:]
#     end
# end
# wadj = MultiStepAdjoint(nw, dis, shot_adj, fidMtx)
# L = dot(vec(shot.p),vec(shot_adj.p)) + dot(vec(shot.vx),vec(shot_adj.vx)) + dot(vec(shot.vz),vec(shot_adj.vz))
# T = dot(w, wadj)

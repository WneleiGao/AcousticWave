using AcousticWave

# generate synthetic data in time domain
function wraperFdModeling()
    nz = 150; nx = 300; ext = 20; iflag = 1;
    dz = 8.; dx = 8.; dt = 1.e-3; f0=20.0; tmax = 3.;
    v  = 3000. * ones(nz, nx);
    # initialize source term
    # first stage
    isz = [40 ,46 ,53 ,54 ,60 ,59 ,67 ,65 ,70 ,69 ,75 ,81 ,80 ,85 ,87 ,92 ,90 ,100,101,110];
    isx = [185,187,190,215,193,212,195,210,200,205,205,196,203,192,209,189,215,186,218,220];
     ot = [0.5,0.4,0.4,0.3,0.3,0.2,0.2,0.1,0.1,0.0,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.5,0.5];
    # second stage
    tmpz = [58 ,61 ,65 ,68 ,72 ,75 ,78 ,82 ,84 ,88 ,92 ,61 ,65 ,68 ,72 ,75 ,78 ,82 ,84 ,88 ];
    tmpx = [140,141,140,139,140,141,139,140,139,140,139,150,151,152,151,150,149,151,150,151];
    tt   = [0.5,0.4,0.4,0.3,0.3,0.2,0.2,0.1,0.1,0.0,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.5,0.5];
    tt   = tt + 1.5
    isz = vcat(isz, tmpz); isx = vcat(isx, tmpx); ot = vcat(ot, tt);

    amp = 1.0*ones(length(isz));
    srcs= InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
    # discretize spatial derivative operator
    vmax = maximum(v); vmin = minimum(v);
    fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, v, f0);
    irz = collect(1:1:nz); irx = 5*ones(Int64, length(irz));
    tmpx = collect(1:1:nx); tmpz = 5*ones(Int64, length(tmpx));
    irz = vcat(irz, tmpz); irx = vcat(irx, tmpx);
    shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax);
    return shot.p
end
@time d = wraperFdModeling();
v = quantile(vec(d),0.95)
imshow(d, cmap="gray", aspect="auto", vmax=v, vmin=-v)

ot = 1.65; nt=1024; dt=0.001; df = 1./(nt*dt);
dc = intercept(d, ot, dt, nt)
imshow(dc, cmap="gray", aspect="auto", vmax=v, vmin=-v)
# # # Source location by stacking FISTA result
# # # ========================================
nz=150; nx=300; ext=50; iflag=1;
dz=0.008; dx=0.008;
M  = Mesh(nt, nz, nx, ext, iflag, dt, dz, dx);
# model parameters
v = 3. * ones(nz, nx); m = 1.0 ./ (v.^2); m = modelExpand(m, M);
# boundary conditions
amp = -1.5*maximum(v)./(ext*dz)*log(0.001); pml= getABC(M, amp);
somm = getSommerfeldBC(M, m);
# laplace operator
L = getLaplacian(M);
# receivers' location
irz  = collect(1:1:nz); irx = 5*ones(Int64, length(irz));
tmpx = collect(1:1:nx); tmpz= 5*ones(Int64, length(tmpx));
irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
# window size and location
woz=10; wox=10; wnz=140; wnx=290; nt=1024;
miniw = 5; maxiw = 55;

rr = setupLocation(miniw, maxiw, M, m, L, pml, somm,
                   irz, irx, woz, wox, wnz, wnx, dc)

loc = Location(M, m, L, pml, somm, irz, irx, woz, wox, wnz, wnx, iw, dc[20,:])
lambda = getlambda(miniw,maxiw,dt,nt,woz,wox,wnz,wnx,irz,irx,M,L,m,pml,somm)


iw = 30; df=0.997; omega = 2*pi*(iw-1)*df
# fdc = fft(dc,1); dobs = fdc[iw, :];
H = getHelmhotzMtx(L, omega, m, pml, somm); H = lufact(H);
A = getAdjointMtx( L, omega, m, pml, somm); A = lufact(A);

# estimate the maximum eigenvalue
lambda = fpower(woz, wox, wnz, wnx, irz, irx, M, H, A, maxiter=60);

lambda = 4.096e-6; mu = 2e-5;
(x40, J) = ffista(dobs, woz, wox, wnz, wnx, irz, irx, M, H, A, mu, lambda, maxiter=100);
tmp40 = reshape(abs(x40), wnz, wnx); imshow(tmp40);

miniw = 10; maxiw= 55;
loc = stackInversion(dc, dt, L, m, pml, somm,
                     woz, wox, wnz, wnx, irz, irx, M, miniw, maxiw)
path = "/Users/wenyue/Desktop/loc1.bin"
(n1,n2,n3) = size(loc);
fid = open(path,"w")
write(fid,n1,n2,n3);
write(fid,vec(loc));
close(fid);

tmp = squeeze(sum(abs(loc),3),3);
v = quantile(vec(tmp), 0.99)
imshow(tmp, vmax=v, vmin=0)



x = collect(1:6);
y = collect(1:6);
rr = Vector{RemoteChannel}(6)
assignment = Vector{Vector{Int64}}(nworkers())
assignment[1] = collect(1:2)
assignment[2] = collect(3:4)
assignment[3] = collect(5:6)
@sync begin
      for (idx, ip) in enumerate(workers())
          @async begin
                 for i = 1 : length(assignment[idx])
                     subidx = assignment[idx][i]
                     rr[subidx] = initRemoteChannel(point, ip, x[subidx], y[subidx])
                 end
          end
      end
end

# # Source location by simply stacking adjoint wavefield
# # ====================================================
# nz=150; nx=300; ext=50; iflag=1;
# dz=0.008; dx=0.008; dt=0.001; nf=2048;
# M  = Mesh(nf, nz, nx, ext, iflag, dt, dz, dx);
# amp=500.; pml= getABC(M, amp);
# m  = 1/9 * ones(nz, nx); m = modelExpand(m, M);
# L  = getLaplacian(M);
# irz  = collect(1:1:nz); irx = 5*ones(Int64, length(irz));
# tmpx = collect(1:1:nx); tmpz= 5*ones(Int64, length(tmpx));
# irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
#
# #shift the recordings
# dc = intercept(d, 0.4, 0.001, 1024)
# SeisPlot(dc, pclip=95, cmap="gray")
#
# miniw=10; maxiw=60; woz=10; wox=10; wnz=140; wnx=290;
# loc = adjoint_stack(dc, dt, miniw, maxiw, L, m, pml, somm, woz, wox, wnz, wnx, irz, irx, M)

# # # Source location by CG method
# # # ====================================================
# nz=150; nx=300; ext=30; iflag=1;
# dz=0.008; dx=0.008; dt=0.001; nf=2048;
# M  = Mesh(nf, nz, nx, ext, iflag, dt, dz, dx);
# amp=500.; pml= getABC(M, amp);
# m  = 1/9 * ones(nz, nx); m = modelExpand(m, M);
# L  = getLaplacian(M);
#
# irz  = collect(1:10:nz); irx = 5*ones(Int64, length(irz));
# tmpx = collect(1:20:nx); tmpz= 5*ones(Int64, length(tmpx));
# irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
#
# #shift the recordings
# dc = intercept(d, 0.300, 0.001, 1200)
# SeisPlot(dc, pclip=95, cmap="gray")
#
# woz=10; wox=10; wnz=140; wnx=290; maxiw=110;
# loc = CG_stack(dc, dt, maxiw, L, m, pml, M,
#                woz, wox, wnz, wnx, irz, irx)
#
# # Source location by IRLS method
# # ====================================================
# nz=150; nx=300; ext=50; iflag=1;
# dz=0.008; dx=0.008; dt=0.001; nf=2048;
# M  = Mesh(nf, nz, nx, ext, iflag, dt, dz, dx);
# amp=500.; pml= getABC(M, amp);
# m  = 1/9 * ones(nz, nx); m = modelExpand(m, M);
# L  = getLaplacian(M);
#
# irz  = collect(1:10:nz); irx = 5*ones(Int64, length(irz));
# tmpx = collect(1:20:nx); tmpz= 5*ones(Int64, length(tmpx));
# irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
#
# #shift the recordings
# dc = intercept(d, 0.300, 0.001, 1200)
# SeisPlot(dc, pclip=95, cmap="gray")
#
# woz=10; wox=10; wnz=140; wnx=290; iw=30;
# df = 1/(dt*size(dc,1)); f = (iw-1)*df; omega= 2*pi*f;
# H = getHelmhotzMtx(L, omega, m, pml); H = lufact(H);
# A = getAdjointMtx( L, omega, m, pml); A = lufact(A);
#
# # dc = fft(dc,1); dobs = dc[iw,:];
# (loc2, J) = CG(dobs, woz, wox, wnz, wnx, irz, irx, M, H, A, Niter=10, mu=0.)
# for i = 1 : 5
# P = goodPass(loc2)
# (loc2, J) = PCG(dobs,loc1,woz,wox,wnz,wnx,irz,irx,M,P,H,A,mu=1e-16)
# end
# imshow(reshape(abs(loc2), wnz, wnx))

# modeling in frequency domai
# =================================================
# using PyPlot, AcousticWave
# nz=150; nx=300; ext=50  ; iflag=1;
# dz=0.008; dx=0.008; dt=0.001; nf=1024;
# M  = Mesh(nf, nz, nx, ext, iflag, dt, dz, dx);
# L  = getLaplacian(M);
# # boundary coefficients
# amp = -1.5*3./(ext*dz)*log(0.001); pml= getABC(M, amp);
# # model parameters
# v = 3. * ones(nz, nx);  m = 1./(v.^2);
# m = modelExpand(m, M);
# somm = getSommerfeldBC(M, m);
# # sources
# isz = [50]; isx = [100]; f0 = 20.; t0=[0.];
# src = FSource(isz, isx, f0, t0, dt, nf);
# # receives
# irz  = collect(1:1:nz); irx = 5*ones(Int64, length(irz));
# tmpx = collect(1:1:nx); tmpz= 5*ones(Int64, length(tmpx));
# irz = vcat(irz, tmpz) ; irx = vcat(irx, tmpx);
# # modeling
# maxiw = 85;
# wfd = getDataTime(src, M, pml, somm, L, maxiw)
# # plot the result
# wfd = real(wfd);
# idx = 330; v = maximum(abs(wfd[idx,:,:]));
# figure(); imshow(wfd[idx,:,:], cmap="seismic", vmax=v, vmin=-v)
# savefig("/Users/wenyue/Desktop/wfd.pdf")
# figure(); plot(wfd[idx,:,150])
# savefig("/Users/wenyue/Desktop/trace.pdf")

# function wraperFdModeling()
#     nz = 100; nx = 200; ext = 20; iflag = 1;
#     dz = 8.; dx = 8.; dt = 1.e-3; f0=20.0; tmax = 1.047;
#     v  = 3000. * ones(nz, nx);
#     # initialize source term
#     isz = [50];
#     isx = [100];
#      ot = [0.0];
#     amp = [1.0];
#     srcs= InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
#     # discretize spatial derivative operator
#     vmax = maximum(v); vmin = minimum(v);
#     fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, v, f0);
#     irz = collect(1:1:nz); irx = 5*ones(Int64, length(irz));
#     tmpx = collect(1:1:nx); tmpz = 5*ones(Int64, length(tmpx));
#     irz = vcat(irz, tmpz); irx = vcat(irx, tmpx);
#     path = "/Users/wenyue/Desktop/wfd.bin"
#     MultiStepForward(path, srcs, fidMtx, tmax=tmax, wtype="p");
# end
# wraperFdModeling()
# spt = readStress(path, idx);
# v = maximum(abs(spt.p));
# imshow(spt.p, cmap="seismic", vmax=v, vmin=-v)

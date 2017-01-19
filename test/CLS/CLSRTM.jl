using PyPlot, AcousticWave

nz = 125; nx = 300; ext = 30; iflag = 1; dz = 20.; dx = 20.; dt = 2e-3; f0 = 10.
v = 2500. * ones(nz, nx); tmax = 4.0;
v[26:100, :] = 3000.
for iz = 1 : 25
    v[50+iz, 151-(iz-1):151+(iz-1)] = 3500.
end
v[101:end, :] = 4000.; vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);

# initialize source term
isx = collect(1:3:nx); ns = length(isx); isz = ones(Int64, ns);
# ot = randST(ns, dt, 0.5);
ot = zeros(ns);
pos_src = hcat(isz, isx); amp = ones(ns);
src = InitSources(pos_src, f0, ot, amp, dt)
# receiver location
irx = collect(1:1:nx); nr = length(irx); irz = ones(Int64, nr);
pos_rec = convert(Array{Int64,2}, hcat(irz, irx))
# position of receiver
shot  = MultiStepForward(pos_rec, src, fidMtx, tmax=tmax);

vda = 2500. * ones(nz, nx);
fidMtx_da = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, vda);
shot_da   = MultiStepForward(pos_rec, src, fidMtx_da, tmax=tmax);
# shot.d[:] = shot.d-shot_da.d
# SeisPlot(shot.d, pclip=98, cmap="PuOr")

# compute source side wavefield with smooth velocity model
v_sm = modSmooth(v, 20)
fidMtx_sm = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v_sm);
pathsrc = "/Users/wenlei/Desktop/src.bin"
MultiStepForward(pathsrc, src, fidMtx_sm, tmax=tmax)

# compute receiver side wavefield with smooth velocity model
pathrec = "/Users/wenlei/Desktop/rec.bin"
MultiStepAdjoint(pathrec, shot, fidMtx_sm)

Ip = RTMimaging_test(pathsrc, pathrec, image_type="p");
SeisPlot(Ip, cmap="PuOr", pclip=99)

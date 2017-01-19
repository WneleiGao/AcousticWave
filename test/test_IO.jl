using PyPlot, AcousticWave
nz = 150; nx = 200; ext= 30  ; iflag = 1 ;
dx = 20.; dz = 20.; dt = 2e-3; tmax  = 4.; f0=10.;
v = 2500.*ones(nz, nx);
v[26:100, :] = 3000.;
for iz = 1 : 25
    v[50+iz, 100-2*(iz-1):100+2*(iz-1)] = 3500.;
end
v[101:125, :] = 4000.; v[126:nz, :] = 4500.;
vl = minimum(v); vu = maximum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vl, vu, v, f0);
# finish testing forward modeling code

isx=collect(30:30:180); isz=ones(Int64,length(isx)); ot=zeros(length(isx)); amp=ones(length(isx));
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
root = homedir(); path = join([root "/Desktop/shot.bin"])
irx = collect(1:1:nx); irz = ones(Int64,length(irx));
shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax);
writeShot!(path, shot)
shotr = readShot(path)
vecnorm(shot.p-shotr.p)

val = quantile(abs(vec(shot.p)), 0.90); imshow(shot.p, vmin=-val, vmax=val, cmap="PuOr", aspect=0.04);

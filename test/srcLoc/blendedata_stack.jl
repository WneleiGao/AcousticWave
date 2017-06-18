using PyPlot, ElasticWave

nz = 100; nx = 200; ext= 20;   iflag = 1;
dx = 20.; dz = 20.; dt = 2e-3; tmax = 4.; f0=10.;
vp = 2500.*ones(nz, nx)   ; vs = 1450.*ones(nz, nx)    ; rho=2.5*ones(nz,nx);
vp[19:68, :     ] = 3000.; vs[19:68 , :     ] = 1752.;
# vp[32:56, 76:125] = 3500.; vs[32:56 , 76:125] = 2020.;
vp[69:85, :     ] = 4000.; vs[69:85 , :     ] = 2309.;
vp[86:nz, :     ] = 4500.; vs[86:nz , :     ] = 2600.;
fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);
root = join([homedir() "/Desktop/shotv"])

flags = vec([false false true true false]);
irx   = collect(1:1:nx); irz = 1*ones(Int64, length(irx));
xs = collect(1:4:nx);
isz = 1; ot = 0.0; isx = 100;
for is = 1 : length(xs)
    isx = xs[is];
    src = InitSource(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);
    shotv = MultiStepForward(irz, irx, src, fidMtx, tmax=tmax);
    path_shotv = join([root "$is" ".bin"])
    WriteShotV(path_shotv, shotv)
end

wfd = ReadWfd(path_wfd, 600)
vx=reshape(wfd.Vx,160,260); vz=reshape(wfd.Vz,160,260);
a=quantile(abs(vec(vx)),0.99); imshow(vx, vmax=a, vmin=-a, cmap="PuOr"); figure();
a=quantile(abs(vec(vz)),0.99); imshow(vz, vmax=a, vmin=-a, cmap="PuOr")

tmp=quantile(abs(vec(shotv.Vx)),0.9); imshow(shotv.Vx, aspect=0.1, vmax=tmp, vmin=-tmp, cmap="PuOr"); figure();
tmp=quantile(abs(vec(shotv.Vz)),0.9); imshow(shotv.Vz, aspect=0.1, vmax=tmp, vmin=-tmp, cmap="PuOr")

root = homedir()
for is  = length(isx)
    src = InitSource(isz, isx[is], nz, nx, ext, iflag, f0, ot, dt, flags);
    irx   = collect(1:1:nx); irz = 1*ones(Int64, length(irx));
    shotv = MultiStepForward(irz, irx, src, fidMtx)
    imshow(shotv.Vx, aspect=0.2, vmax=maximum(abs(shotv.Vx)), vmin=-maximum(abs(shotv.Vx)), cmap="PuOr"); figure();
    imshow(shotv.Vz, aspect=0.2, vmax=maximum(abs(shotv.Vz)), vmin=-maximum(abs(shotv.Vz)), cmap="PuOr")
    removeDirect(shotv, vp[1], src, dx, dt)
    figure(); imshow(shotv.Vx, aspect=0.2, vmax=maximum(abs(shotv.Vx)), vmin=-maximum(abs(shotv.Vx)), cmap="PuOr");
    figure(); imshow(shotv.Vz, aspect=0.2, vmax=maximum(abs(shotv.Vz)), vmin=-maximum(abs(shotv.Vz)), cmap="PuOr");
    path_shotv = join([root "/Desktop/shot" "$is"])
    WriteShotV(path_shotv)
end

function removeFirstArrival(isx)
    nz = 30; nx = 203; ext= 30;   iflag = 1;
    dx = 5. ; dz = 5. ; dt = 4e-4; tmax = 0.5; f0=30.;
    vp =3500.*ones(nz,nx);  vs=1750.*ones(nz,nx); rho=2.5*ones(nz,nx);
    fidMtx = CreateFidMtx(nz, nx, ext, iflag, vp, vs, rho, dz, dx, dt, f0);

    isz = 1; ot = zeros(length(isx));
    flags = vec([false false true true false]);
    srcs  = InitMultiSources(isz, isx, nz, nx, ext, iflag, f0, ot, dt, flags);

    irx   = collect(1:1:nx); irz = 1*ones(Int64, length(irx));
    shotv = MultiStepForward(irz, irx, srcs, fidMtx)
    return shotv
end

imshow(shotv.Vx, aspect=0.1); figure(); imshow(shotv.Vz, aspect=0.1)
root = homedir()
path_shotv = join([root "/Desktop/shotv.bin"])
WriteShotV(path_shotv, shotv);

shotv_FA = removeFirstArrival();
shotv.Vx = shotv.Vx - shotv_FA.Vx;
shotv.Vz = shotv.Vz - shotv_FA.Vz;

path_shotv_DFA = join([root "/Desktop/shotv_DFA.bin"])
WriteShotV(path_shotv_DFA, shotv);

path_pv = join([root "/Desktop/pv.bin"])
MultiStepForward(path_pv, srcs, fidMtx, dz, dx)

(dm, du) = MultiStepAdjoint(shotv, fidMtx, path_pv);
imshow(dm); figure(); imshow(du);

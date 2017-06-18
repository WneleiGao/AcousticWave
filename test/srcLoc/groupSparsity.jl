using PyPlot, AcousticWave
# =============the result of adjoint================================
nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 8e-4;
v  = 3000. * ones(nz, nx);
# initialize source term
pos_src = [101 75; 101 125]; f0 = 30.0; ot = zeros(size(pos_src,1));
src = InitSources(pos_src, f0, ot, dt);
(tl, tu) = SrcRange(src)
ntw = round(Int64, (tu-tl)/dt+1)

vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
iz = collect(5:5:201); ix = 20 * ones(Int64, length(iz));
pos_rec = hcat(iz, ix);
shot = MultiStepForward(pos_rec, src, fidMtx, tmax=0.5);

path_x = "/Users/wenlei/Desktop/x.bin"
path_s = "/Users/wenlei/Desktop/s.bin"
path_p = "/Users/wenlei/Desktop/p.bin"
W = CGjoint_GP(ntw, path_x, path_s, path_p, shot, fidMtx, pos_rec, mu=10.0, outer=3)
SeisPlot(readSrcSpt(path_x, 10))

# check data fitting
shot_pre = MultiStepForward(pos_rec, path_x, fidMtx, tmax=0.5)
it = 40
plot(shot.d[:,it])
plot(shot_pre.d[:,it])

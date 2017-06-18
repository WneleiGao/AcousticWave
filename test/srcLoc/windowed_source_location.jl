using PyPlot, AcousticWave

nz = 122; nx = 384; ext = 20; iflag = 1;
dz = 10.; dx = 10.; dt = 1.0e-3; f0=12.0; tmax = 1.5;
root = homedir(); path = join([root "/Desktop/marm.dat"])
v = readdlm(path); v = reshape(v,nz,nx);
v = modSmooth(v,10); vl=minimum(v); vu=maximum(v);
SeisPlot(v, hbox=4.5, wbox=12, dx=10, dy=10, xlabel="x (m)", ylabel="z (m)", cmap="jet", vmax=vu, vmin=vl)
colorbar(); savefig("/Users/wenyue/Desktop/vel.pdf");
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vu, vl, v, f0);

isz = [81, 81, 81, 81, 81]; isx = [140, 170, 200, 230, 260];
ot  = [0.0, 0.04, 0.08, 0.10, 0.13]; amp = [1.0, -1.0, 1.0, -1.0, 1.0];
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);

irx  = collect(1:3:nx); irz = 1 * ones(Int64, length(irx));
shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax)
oz = 1; ox = 1; wnz = nz; wnx = nx; wnt=500;
wsc_adj = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, shot, fidMtx);
tmp = zeros(wnt, 5);
tmp[:,1] = vec(wsc_adj.p[81,140,:])
tmp[:,2] = vec(wsc_adj.p[81,170,:])
tmp[:,3] = vec(wsc_adj.p[81,200,:])
tmp[:,4] = vec(wsc_adj.p[81,230,:])*0.5
tmp[:,5] = vec(wsc_adj.p[81,260,:])*0.4
# l = maximum(tmp)
# for ix = 1 : size(tmp,2)
#     u = maximum(abs(tmp[:,ix]))
#     r = u / l
#     tmp[:,ix] = tmp[:,ix] / r
# end
SeisPlot(tmp, style="wiggles",wbox=5, hbox=5, xcur=0.8, ox=1, dx=1, dy=dt, xlabel="Trace number", ylabel="Time (s)", name="/Users/wenyue/Desktop/adjoint.pdf")
wd = engdis(wsc_adj); whd = zeros(nz,nx);
whd[oz:oz+wnz-1, ox:ox+wnx-1] = wd[:,:]
SeisPlot(whd, wbox=12, hbox=5, dx=dx, dy=dz, xlabel="x (m)", ylabel="z (m)", name="/Users/wenyue/Desktop/dis_adjont.pdf")


oz = 61; ox = 111; wnz = 40; wnx = 180;
wsc = Srcs2Wsc(oz, ox, wnz, wnx, srcs); wnt = 500;
wd = engdis(wsc);
wd[81-2:81+2, 140-2:140+1] = 1.
wd[81-2:81+2, 170-2:170+1] = 1.
wd[81-2:81+2, 200-2:200+1] = 1.
wd[81-2:81+2, 230-2:230+1] = 1.
wd[81-2:81+2, 260-2:260+1] = 1.
SeisPlot(wd, wbox=12, hbox=5, dx=dx, dy=dz, xlabel="x (m)", ylabel="z (m)", name="/Users/wenyue/Desktop/dis_true.pdf")


irx  = collect(1:3:nx); irz = 1 * ones(Int64, length(irx));
shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax)
SeisPlot(shot.p, pclip=99.7, dx=1, dy=dt, xlabel="Trace number", ylabel="Time (s)", name="/Users/wenyue/Desktop/shot.pdf")

mu = 0.48; lambda = 380.0
(J, wsc2) = FISTA(oz, ox, wnz, wnx, wnt, shot, fidMtx, mu, lambda, niter=30)
wd = engdis(wsc2); whd = zeros(nz,nx);
whd[oz:oz+wnz-1, ox:ox+wnx-1] = wd[:,:]
SeisPlot(whd, wbox=12, hbox=5, dx=dx, dy=dz, xlabel="Length (m)", ylabel="Depth (m)", name="/Users/wenyue/Desktop/dis_inv.pdf")

wsc2 = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, shot, fidMtx)
tmp = zeros(wnt, 5);
tmp[:,1] = vec(wsc2.p[21,30,:])
tmp[:,2] = vec(wsc2.p[21,60,:])
tmp[:,3] = vec(wsc2.p[21,90,:])
tmp[:,4] = vec(wsc2.p[21,120,:])
tmp[:,5] = vec(wsc2.p[21,150,:])
l = maximum(tmp)
for ix = 1 : size(tmp,2)
    u = maximum(abs(tmp[:,ix]))
    r = u / l
    tmp[:,ix] = tmp[:,ix] / r
end
SeisPlot(tmp, style="wiggles",wbox=5, hbox=5, xcur=0.8, ox=1, dx=1, dy=dt, xlabel="Trace number", ylabel="Time (s)", name="/Users/wenyue/Desktop/adjoint.pdf")
wd = engdis(wsc2); whd = zeros(nz,nx);
whd[oz:oz+wnz-1, ox:ox+wnx-1] = wd[:,:]
SeisPlot(whd, wbox=12, hbox=5, dx=dx, dy=dz, xlabel="x (m)", ylabel="z (m)", name="/Users/wenyue/Desktop/dis_true.pdf")


tmp1 = zeros(500, 5); wnt = wsc.nt
tmp1[1:wnt,1] = vec(wsc.p[21,30,:])
tmp1[1:wnt,2] = vec(wsc.p[21,60,:])
tmp1[1:wnt,3] = vec(wsc.p[21,90,:])
tmp1[1:wnt,4] = vec(wsc.p[21,120,:])
tmp1[1:wnt,5] = vec(wsc.p[21,150,:])
SeisPlot(tmp1, style="wiggles", xcur=0.8, ox=1, dx=1, dy=dt, xlabel="Trace number", ylabel="Time (s)", name="/Users/wenyue/Desktop/true.pdf")


# =======compute maximum eigenvalue=========
nt = shot.nt
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
SeisPlot(shot1.p)
wsc = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, shot1, fidMtx);

wsc = RandWsc(oz, ox, wnz, wnx, ext, iflag, dt, wnt, f0)
lambda = power_wsc(irz, irx, fidMtx, wsc, tmax, niter=40)

function RandWsc(oz::Int64, ox::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, nt::Int64, f0::Float64)
    w = Ricker(f0, dt)
    l = floor(Int64, length(w)/2)
    a = randn(nz, nx, nt)
    for iz = 1 : nz
        for ix = 1 : nx
            tmp = vec(a[iz,ix,:])
            tmp1= conv(tmp, w)
            a[iz,ix,:] = tmp1[l+1:end-l]
        end
    end
    wsc = WSrcCube(oz, ox, nz, nx, ext, iflag, 0.0, dt, nt, a)
    return wsc
end

# pass the forward modelling test and dot product test.
function power_wsc(irz::Array{Int64,1}, irx::Array{Int64,1}, fidMtx::FidMtx, wsc::WSrcCube, tmax::Float64; niter=20)
    # compute the dominant eigenvalue of A'A
    oz = wsc.oz; ox = wsc.ox; wnz = wsc.nz; wnx = wsc.nx; wnt =wsc.nt;
    lambda = 0.0
    for iter = 1 :niter
        shot = MultiStepForward(irz, irx, wsc, fidMtx, tmax=tmax)
        # shot.p = shot.p / sqrt(680)
        wsc  = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, shot, fidMtx)
        # wsc.p = wsc.p / sqrt(680)
        lambda = vecnorm(wsc.p)
        wsc.p  = wsc.p / lambda
        println("iteration: $iter, maximum eig: $lambda")
    end
    return lambda
end

function engdis(wsc::WSrcCube)
    (nz, nx, nt) = size(wsc.p)
    w = zeros(nz, nx)
    for ix = 1 : nx
        for iz = 1 : nz
            for it = 1 : nt
                w[iz, ix] = w[iz,ix] + wsc.p[iz,ix,it] * wsc.p[iz,ix,it]
            end
            w[iz,ix] = sqrt(w[iz,ix])
        end
    end
    return w
end

function softThresh_GP!(wsc::WSrcCube, T::Float64)
    (nz, nx, nt) = size(wsc.p)
    w = zeros(nz, nx)
    for ix = 1 : nx
        for iz = 1 : nz
            for it = 1 : nt
                w[iz, ix] = w[iz,ix] + wsc.p[iz,ix,it] * wsc.p[iz,ix,it]
            end
            w[iz,ix] = sqrt(w[iz,ix])
            if w[iz,ix] < T
               w[iz,ix] = 0.
            else
               w[iz,ix] = (w[iz,ix]-T) / w[iz,ix]
            end
        end
    end
    for it = 1 : nt
        wsc.p[:,:,it] = wsc.p[:,:,it] .* w
    end
    return nothing
end

function softThresh_L1!(wsc::WSrcCube, T::Float64)
    (nz, nx, nt) = size(wsc.p)
    for it = 1 : nt
        for ix = 1 : nx
            for iz = 1 : nz
                if abs(wsc.p[iz,ix,it]) < T
                   wsc.p[iz,ix,it] = 0
                else
                   wsc.p[iz,ix,it] = (abs(wsc.p[iz,ix,it])-T)*sign(wsc.p[iz,ix,it])
                end
            end
        end
    end
    return nothing
end

function FISTA(oz::Int64, ox::Int64, wnz::Int64, wnx::Int64, wnt::Int64, d::Shot, fidMtx::FidMtx, mu::Float64, lambda::Float64; niter=20)
    J = zeros(niter + 1)
    T = mu / (2*lambda)
    x0 = InitWSrcCube(oz, ox, wnz, wnx, d.ext, d.iflag, 0.0, d.dt, wnt)
    dtmp = InitShot(0, 0, d.nz, d.nx, d.ext, d.iflag, d.irz, d.irx, 0.0, d.dt, d.nt)
    xk = CopyWsc(x0)
    yk = CopyWsc(x0)
    t0 = 1.0
    irz = d.irz; irx = d.irx; tmax=(d.nt-1)*d.dt
    for k = 1 : niter
        Gy = MultiStepForward(irz, irx, yk, fidMtx, tmax=tmax)
        dtmp.p = shot.p - Gy.p
        cost  = (vecnorm(dtmp.p))^2
        println("iterations: $k, cost $cost"); J[k] = cost;
        mtmp = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, dtmp, fidMtx)
        xk.p = yk.p + mtmp.p/lambda
        softThresh_GP!(xk, T)
        t = (1 + sqrt(1+4*t0^2)) / 2
        yk.p = xk.p + (t0-1)/t*(xk.p-x0.p)
        t0 = t
        x0.p = copy(xk.p)
    end
    return J, xk
end


# =================simple model ============================
using PyPlot, AcousticWave
nz = 201; nx = 201; ext = 20; iflag = 1;
dz = 15.; dx = 15.; dt = 1e-3; f0=10.0; tmax = 1.5;
v  = 3000. * ones(nz, nx); vl=minimum(v); vu=maximum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vu, vl, v, f0);

isz = [75,125]; isx = [101,101];  ot = [0.0,0.0]; amp = [1.0,-1.0]
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);

oz = 45; ox = 71; wnz = 110; wnx = 60;
wsc = Srcs2Wsc(oz, ox, wnz, wnx, srcs); wnt = wsc.nt;

irz = collect(10:5:190); irx = 20 * ones(Int64, length(irz));
shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax, otype="p");
SeisPlot(shot.p, pclip=99.7, cmap="PuOr", dy=dt, dx=1, xlabel="Trace Number", ylabel="Time (s)", name="/Users/wenyue/Desktop/shot1.pdf")

oz = 1; ox = 1; wnz = 200; wnx = 200;
wsc1 = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, shot, fidMtx)
w = engdis(wsc1)
dis = zeros(nz, nx)
dis[oz:oz+wnz-1, ox:ox+wnx-1] = w[:,:]

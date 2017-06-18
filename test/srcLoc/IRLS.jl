using PyPlot, AcousticWave

nz = 122; nx = 384; ext = 20; iflag = 1;
dz = 10.; dx = 10.; dt = 1.0e-3; f0=12.0; tmax = 1.5;
root = homedir(); path = join([root "/Desktop/marm.dat"])
v = readdlm(path); v = reshape(v,nz,nx);
v = modSmooth(v,10); vl=minimum(v); vu=maximum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vu, vl, v, f0);

isz = [81, 81, 81, 81, 81]; isx = [140, 170, 200, 230, 260];
ot  = [0.0, 0.04, 0.08, 0.10, 0.13]; amp = [1.0, 1.0, 1.0, 1.0, 1.0];
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);

oz = 51; ox = 110; wnz = 60; wnx = 190;
wsc = Srcs2Wsc(oz, ox, wnz, wnx, srcs); wnt = 500;

irx  = collect(1:3:nx); irz = 1 * ones(Int64, length(irx));
shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax);

wsc1 = CGLS(oz, ox, wnz, wnx, wnt, shot, fidMtx)

function CGLS(oz, ox, wnz, wnx, wnt, d, fidMtx; mu=0.1, niter=10)
    tmax = (d.nt-1) * d.dt
    cost = Float64[]
    r = InitShot(0, 0, d.nz, d.nx, d.ext, d.iflag, d.irz, d.irx, d.ot,  d.dt, d.nt)
    r.p = copy(d.p)
    g = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, r, fidMtx)
    m = InitWSrcCube(oz, ox, wnz, wnx, d.ext, d.iflag, d.ot, d.dt, wnt)
    s = InitWSrcCube(oz, ox, wnz, wnx, d.ext, d.iflag, d.ot, d.dt, wnt); s.p[:] = g.p[:];
    gamma = l2Wsc(g)
    gamma0= gamma
    cost0 = (vecnorm(r.p))^2
    push!(cost, 1.0)
    for iter = 1 : niter
        t = MultiStepForward(d.irz, d.irx, s, fidMtx, tmax=tmax)
        delta = (vecnorm(t.p))^2 + mu*l2Wsc(s)
        alpha = gamma / delta
        m.p   = m.p + alpha*s.p
        r.p   = r.p - alpha*t.p
        g     = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, r, fidMtx)
        g.p   = g.p - mu*m.p
        gamma0= copy(gamma)
        gamma = l2Wsc(g)
        cost1 = (vecnorm(r.p))^2 + mu*l2Wsc(m)
        tmp = cost1/cost0
        println("iter $iter, cost $tmp")
        push!(cost, cost1/cost0)
        beta = gamma / gamma0
        s.p = g.p + beta*s.p
    end
    return m, cost
end

function l2Wsc(wsc::WSrcCube)
    (nz, nx, nt) = size(wsc.p)
    tmp = 0.0
    for it =1 : nt
        for ix = 1 : nx
            for iz = 1 : nz
                tmp += wsc.p[iz,ix,it]*wsc.p[iz,ix,it]
            end
        end
    end
    return tmp
end









n = 4;
sg = divideGroup(shot, n)
wscg = Array(WSrcCube, n)
for i = 1 : n
    wscg[i] = MultiStepAdjoint(oz, ox, wnz, wnx, wnt, sg[i], fidMtx)
end
I = ones(wnz, wnx, wnt);
for it = 1 : wnt
    for ig = 1 : n
        I[:,:,it] = I[:,:,it] .* wscg[ig].p[:,:,it]
    end
end
Ic = zeros(wnz, wnx)
for it = 1 : wnt
    Ic += I[:,:,it] .* I[:,:,it]
end

function divideGroup(shot::Shot, n::Int64)
    gn = floor(Int64, shot.nr/n)
    R  = Array(Shot, n)
    for i = 1 : n-1
        irz = shot.irz[(i-1)*gn+1:i*gn]
        irx = shot.irx[(i-1)*gn+1:i*gn]
        r = InitShot(0, 0, shot.nz, shot.nx, shot.ext, shot.iflag, irz, irx, shot.ot, shot.dt, shot.nt)
        r.p[:,1:gn] = shot.p[:,(i-1)*gn+1:i*gn]
        R[i] = r
    end
    irz = shot.irz[(n-1)*gn+1:end]
    irx = shot.irx[(n-1)*gn+1:end]
    r = InitShot(0, 0, shot.nz, shot.nx, shot.ext, shot.iflag, irz, irx, shot.ot, shot.dt, shot.nt)
    r.p[:,1:end] = shot.p[:,(n-1)*gn+1:end]
    R[n] = r
    return R
end

function randGroup(shot::Shot, n::Int64)
    ind = randperm(shot.nr)
    gn = floor(Int64, shot.nr/n)
    R  = Array(Shot, n)
    for ig = 1 : n
        if ig <= n-1
           il = (ig-1)*gn+1
           iu = ig*gn
        elseif ig == n
           il = (n-1)*gn+1
           iu = shot.nr
        end
        p  = zeros(shot.nt, iu-il+1)
        irz = zeros(Int64, iu-il+1)
        irx = zeros(Int64, iu-il+1)
        for it = il : iu
            irz[it-il+1] = shot.irz[ind[it]]
            irx[it-il+1] = shot.irx[ind[it]]
            p[:,it-il+1] = shot.p[:,ind[it]]
        end
        r = InitShot(0, 0, shot.nz, shot.nx, shot.ext, shot.iflag, irz, irx, shot.ot, shot.dt, shot.nt)
        r.p[:,:] = p[:,:]
        R[ig] = r
    end
    return R
end

function slidenormalize(I::Array{Float64,3}, f0, dt)
    w = Ricker(f0, dt)
    hw = floor(Int64, length(w)/2)
    (nz, nx, nt) = size(I)
    for ix = 1 : nx
        for iz = 1 : nz
            for it = 1 : nt
                ql = 1  > it-hw ? 1  : it-hw
                qu = nt < it+hw ? nt : it+hw
                I[iz,ix,it] = I[iz,ix,it] / (maximum(abs(vec(I[iz,ix,ql:qu])))+1e-10)
            end
        end
    end
    return nothing
end

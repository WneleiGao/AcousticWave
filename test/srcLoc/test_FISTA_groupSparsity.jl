# estimate largest eigenvalue based on power's method
using PyPlot, AcousticWave
# =============the result of adjoint================================
nz = 151; nx = 151; ext = 20; iflag = 1;
dz = 20.; dx = 20.; dt = 1e-3; f0 = 14.0;
v  = 3000. * ones(nz, nx); tmax = 1.5;

# set up mmodel
vu = maximum(v); vl = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vu, vl, v, f0);

# initialize source term
isz = [75,125]; isx = [101,101];  ot = [0.0,0.0]; amp = [1.0,-1.0]
srcs = InitMultiSources(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
(tl, tu) = SrcRange(srcs)
ntw = round(Int64, (tu-tl)/dt+1)

# receivers location
irz = collect(10:10:nz); irx = 20 * ones(Int64, length(irz));
shot = MultiStepForward(irz, irx, srcs, fidMtx, tmax=tmax);

# power's method to determine maximum eigenvalue
nt = shot.nt;
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
root = homedir(); path = join([root "/Desktop/power.bin"]);
MultiStepAdjoint(path, ntw, shot1, fidMtx)
lambda = power_wave(irz, irx, ntw, path, tmax, niter=20);
# for this setting lambda = 3371.74


path_x = "/Users/wenlei/Desktop/x.bin"
path_x0 = "/Users/wenlei/Desktop/x0.bin"
path_y = "/Users/wenlei/Desktop/y.bin"
path_tmp = "/Users/wenlei/Desktop/tmp.bin"
mu = 0.2; lambda = 380.; max_iter = 20; tmax=0.5;
FISTA(ntw, path_x, path_x0, path_y, path_tmp, shot, fidMtx, tmax, mu, lambda, max_iter)
SeisPlot(readSrcSpt(path_x, 10))

function power_wave(irz::Array{Int64,1}, irx::Array{Int64,1}, ntw::Int64, path_m::ASCIIString, tmax::Float64; niter=20)
    lambda = 0.0
    for iter = 1 :niter
        d_shot = MultiStepForward(irz, irx, path_m, fidMtx, tmax=tmax)
        MultiStepAdjoint(path_m, ntw, d_shot, fidMtx)
        lambda = L2normStress(path_m)
        scaleStress(path_m, 1/lambda)
        println("iteration: $iter, maximum eig: $lambda")
    end
    return lambda
end


function softThresh_GP!(path_x, path_y::ASCIIString, T::Float64)
    eng = engSrcSpt(path_y)
    mask = zeros(eng)
    for i = 1 : length(eng)
        mask[i] = eng[i] - T
        if mask[i] < 0.0
           mask[i] = 0.0
        else
           mask[i] = mask[i] / eng[i]
        end
    end
    (nz, nx, ot, dt, nt) = InfoSrcSpt(path_y)
    fid_x = open(path_x, "w"); fid_y = open(path_y, "r");
    write(fid_x, nz, nx, ot, dt, nt)
    pre = sizeof(Float64)*5
    lspt = sizeof(Float64)*nz*nx
    for it = 1 : nt
        position = pre + (it-1)*lspt
        seek(fid_y, position)
        m = read(fid_y, Float64, nz*nx)
        m = m .* mask
        seek(fid_x, position)
        write(fid_x, m)
    end
    close(fid_y)
    close(fid_x)
    return nothing
end

function updatey!(path_y::ASCIIString, path_x::ASCIIString, path_x0::ASCIIString, t0::Float64, t::Float64)
    (nz, nx, ot, dt, nt) = InfoSrcSpt(path_x)
    fid_y = open(path_y, "w"); fid_x = open(path_x, "r"); fid_x0 = open(path_x0, "r");
    write(fid_y, nz, nx, ot, dt, nt)
    pre = sizeof(Float64)*5
    lspt= sizeof(Float64)*nz*nx
    c = (t0-1.) / t
    for it = 1 : nt
        position = pre + (it-1)*lspt
        seek(fid_x, position); seek(fid_x0, position); seek(fid_y, position);
        x = read(fid_x, Float64, nz*nx)
        x0= read(fid_x0,Float64, nz*nx)
        y = x + c*(x-x0)
        write(fid_y, y)
    end
    close(fid_y); close(fid_x); close(fid_x0);
    return nothing
end

function FISTA(ntw::Int64, path_x::ASCIIString, path_x0::ASCIIString, path_y::ASCIIString, path_tmp::ASCIIString, shot::Shot, fidMtx::FidMtx, tmax::Float64, mu::Float64, lambda::Float64, max_iter::Int64)
    J = zeros(max_iter+1)
    T = mu / (2*lambda)
    cost = dot(vec(shot.d), vec(shot.d)); J[1] = cost;
    println("iterations: 0, object value $cost")
    MultiStepAdjoint(path_x0, ntw, shot, fidMtx)
    cp(path_x0, path_y, remove_destination=true)
    t0 = 1.0
    irz = shot.irz; irx = shot.irx;
    for k = 1 : max_iter
        shot_tmp = MultiStepForward(irz, irx, path_y, fidMtx, tmax=tmax)
        shot_tmp.d = shot.d - shot_tmp.d
        cost = dot(vec(shot_tmp.d), vec(shot_tmp.d)) + mu*sum(engSrcSpt(path_y))
        println("iterations: $k, object value $cost"); J[k+1] = cost;
        MultiStepAdjoint(path_tmp, ntw, shot_tmp, fidMtx)
        AddSrcSpt(path_y, path_tmp, 1/lambda)
        softThresh_GP!(path_x, path_y, T)
        t = (1 + sqrt(1+4*t0^2)) / 2
        updatey!(path_y, path_x, path_x0, t0, t)
        t0 = t; cp(path_x, path_x0, remove_destination=true)
    end
    return J
end

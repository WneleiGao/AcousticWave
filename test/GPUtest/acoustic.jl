using PyPlot, AcousticWave

nz = collect(200:50:600); nx = collect(200:50:600);

function testCpuTime(nz, nx)
    ext = 20; iflag = 1;
    dz = 15.; dx = 15.; dt = 1.e-3; f0=10.0; tmax = 2.0;
    v  = 2500. * ones(nz, nx);

    vmax = maximum(v); vmin = minimum(v);
    fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, v, f0);

    # initialize source term
    isz = 5; isx = 100; ot = 0.; amp = 1.0;
    src = InitSource(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);

    irx = collect(1:2:nx); irz = 5*ones(Int64, length(irx));

    tic()
    shot = MultiStepForward(irz, irx, src, fidMtx, tmax=tmax);
    tc = toq()
    return tc
end

function testGpuTime(nz, nx)
    ext = 20; iflag = 1;
    dz = 15.; dx = 15.; dt = 1.e-3; f0=10.0; tmax = 2.0;
    v  = 2500. * ones(nz, nx);

    vmax = maximum(v); vmin = minimum(v);
    fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, v, f0);

    # initialize source term
    isz = 5; isx = 100; ot = 0.; amp = 1.0;
    src = InitSource(isz, isx, nz, nx, ext, iflag, ot, dt, f0, amp);
    irx = collect(1:2:nx); irz = 5*ones(Int64, length(irx));

    gm = SendFidMtx2GPU(fidMtx);
    tic()
    shot = MultiStepForwardGPU(irz, irx, src, gm, tmax=tmax);
    tg = toq()
    return tg
end

n = length(nz)
cputime = zeros(n)
gputime = zeros(n)
for i = 1 : n
    println("$i")
    cputime[i] = testCpuTime(nz[i], nx[i])
    gputime[i] = testGpuTime(nz[i], nx[i])
end

plot(nz, cputime, linewidth=2, label="cpuTime")
plot(nz, gputime, linewidth=2, label="gpuTime")
legend(fontsize=15, loc = "upper left" )
ax = gca()
plt[:setp](ax[:get_xticklabels](), fontsize=15)
plt[:setp](ax[:get_yticklabels](), fontsize=15)
xlabel("Size of Model", fontsize=15)
ylabel("Time (s)", fontsize=15)



plot(nz, cputime ./ gputime, linewidth=2)
ax = gca()
plt[:setp](ax[:get_xticklabels](), fontsize=15)
plt[:setp](ax[:get_yticklabels](), fontsize=15)
xlabel("Size of Model", fontsize=15)
ylabel("Speed up", fontsize=15)

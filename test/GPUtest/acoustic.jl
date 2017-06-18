using PyPlot, AcousticWave

nz = collect(200:50:600); nx = collect(200:50:600);
cputime = zeros(length(nz)); gputime = zeros(length(nz))

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

    gm = SendFidMtx2GPU(fidMtx);
    tic()
    shot = MultiStepForwardGPU(irz, irx, src, gm, tmax=tmax);
    tg = toq()
    return tg
end

for i = 1 : length(n)


end

function waveAnim(path::String , pathout="NULL"; inter=5, wbox=6, hbox=6, ox=0.0, dx=5.0, oz=0.0, dz=5.0, pclip=100, cmap="PuOr", aspect="auto", interval=40)
    fid = open(path, "r")
    nz  = read(fid, Int64); nx    = read(fid, Int64);
    ext = read(fid, Int64); iflag = read(fid, Int64);
    dt  = read(fid, Float64);
    head    = sizeof(Int64) * 5
    sptSize = sizeof(Float64) * nz * nx
    nt =  round(Int64, (filesize(fid)-head)/sptSize)
    fig = plt.figure(figsize=(wbox, hbox))
    ims = PyCall.PyObject[]
    vmax = 0.2; vmin=-0.2;
    for it = 1 : inter : nt
        position = head + (it-1)*sptSize; seek(fid, position);
        d = reshape(read(fid, Float64, nz*nx), nz, nx)
        # vmax = quantile(vec(abs(d)), pclip/100.)
        # vmin = -vmax
        im = plt.imshow(d, cmap=cmap, vmin=vmin,vmax=vmax,extent=[ox, ox+(size(d,2)-1)*dx, oz+(size(d,1)-1)*dz,oz], aspect=aspect)
        push!(ims, PyCall.PyObject[im])
    end
    close(fid)
    ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
    pathout = join([pathout ".mp4"])
    ani[:save](pathout)
    return nothing
end

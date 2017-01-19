function RTMimage(I::Array{Float64,2}, spt::SnapShot, path::String )
    nz = spt.nz; nx = spt.nx; ext = spt.ext; iflag = spt.iflag;
    it = spt.it
    if spt.iflag == 1
       Nz = nz + 2*ext
       upper = ext
    elseif spt.iflag == 2
       Nz = nz +   ext
       upper = 0
    end
    Nx = nx + 2*ext
    sts= readStress(path, it)
    p  = reshape(spt.pz+spt.px, Nz, Nx) * 1/2
    p  = p[upper+1:upper+nz, ext+1:ext+nx]
    I  = I + p .* sts.p
    return I
end

function RTMimaging(pathsrc::String , pathrec::String ; image_type="p")
    (nz, nx, ext, iflag, dt, nt) = InfoStress(pathsrc)
    fidsrc = open(pathsrc, "r");
    fidrec = open(pathrec, "r");
    head   = sizeof(Int64)*5
    if image_type == "p"
       Ip = zeros(nz, nx)
       lsts = size(Int64)*nz*nx
       for it = 1 : nt
           position = head + (it-1)*lsts
           seek(fidsrc, position); seek(fidrec, position)
           ps = reshape(read(fidsrc, Float64, nz*nx), nz, nx)
           pr = reshape(read(fidrec, Float64, nz*nx), nz, nx)
           Ip = Ip + ps .* pr
       end
       close(fidsrc); close(fidrec);
       return Ip
    elseif image_type == "vp"
       Ivz = zeros(nz, nx); Ivx = zeros(nz, nx); Ip = zeros(nz, nx);
       lwfd = size(Int64)*nz*nx*3
       for it = 1 : nt
           position = head + (it-1)*lwfd
           seek(fidsrc, position); seek(fidrec, position)
           Svz = reshape(read(fidsrc, Float64, nz*nx), nz, nx)
           Svx = reshape(read(fidsrc, Float64, nz*nx), nz, nx)
           Sp  = reshape(read(fidsrc, Float64, nz*nx), nz, nx)
           Rvz = reshape(read(fidrec, Float64, nz*nx), nz, nx)
           Rvx = reshape(read(fidrec, Float64, nz*nx), nz, nx)
           Rp  = reshape(read(fidrec, Float64, nz*nx), nz, nx)
           Ivz = Ivz + Svz .* Rvz
           Ivx = Ivx + Svx .* Rvx
           Ipz = Ipz + Sp  .* Rp
       end
       close(fidsrc); close(fidrec);
       return Ivz, Ivx, Ip
    end
end

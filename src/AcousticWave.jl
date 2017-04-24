module AcousticWave
    using Requires
    const spmatveclib = abspath(joinpath(splitdir(Base.source_path())[1],"..","deps","builds","spmatvec"))
    include("dataStructure/dataStructure.jl")
    include("frequencyWave/frequencyWave.jl")
    include("plot/plotting.jl")
    include("imaging/imaging.jl")
    include("solver/solver.jl")
    include("spmatvec/spmatvec.jl")
    @require PyPlot include("plot/plotting.jl")
end

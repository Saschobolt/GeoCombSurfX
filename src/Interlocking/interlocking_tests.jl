using Logging
include("../Polyhedron.jl")
include("interlocking.jl")

io = open("log.txt", "w+")
logger = SimpleLogger(io)
global_logger(logger)
for filename in ["[-05,1].jl", "[-05,05].jl", "[-07,08].jl", "[-025,075].jl", "[05,0].jl", "[05,05].jl", "[06,06].jl", "[025,0].jl", "[025,025].jl"]
    include(filename)
    write(io, "\n $(filename) \n")
    ass = [Polyhedron(facets = faces[i], verts = coordinates[i]) for i in eachindex(faces)]
    frame = collect(2:9)

    write(io, "wangtest \n")
    try
        model = wangtest(ass, frame)
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nassemblability \n")
    write(io, "Blocks 2,3 rausnehmen\n")
    inds = collect(1:9)
    assembly = ass[inds]
    fr = Int[1,4,5,6,7,8,9]
    try
        model = wangtest(assembly, fr)
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nBlocks 8,9 rausnehmen\n")
    inds = [1,4,5,6,7,8,9]
    assembly = ass[inds]
    fr = Int[1,2,3,4,5]
    try
        model = wangtest(assembly, fr)
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nBlocks 1,7 rausnehmen\n")
    inds = [1,4,5,6,7]
    assembly = ass[inds]
    fr = Int[2,3,4]
    try
        model = wangtest(assembly, fr)
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nBlocks 4 rausnehmen\n")
    inds = [4,5,6]
    assembly = ass[inds]
    fr = Int[2,3]
    try
        model = wangtest(assembly, fr)
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nBlocks 5 rausnehmen\n")
    inds = [5,6]
    assembly = ass[inds]
    fr = Int[1]
    try
        model = wangtest(assembly, fr)
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "wangtest only translational motions\n")
    try
        model = wangtest(ass, frame)
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nassemblability in plane\n")
    write(io, "Blocks 2,3 rausnehmen\n")
    inds = collect(1:9)
    assembly = ass[inds]
    fr = Int[1,4,5,6,7,8,9]
    try
        model = wangtest(assembly, fr; basis_translation = [1 0; 0 1; 0 0], basis_rotation = zeros(3,2))
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nBlocks 8,9 rausnehmen\n")
    inds = [1,4,5,6,7,8,9]
    assembly = ass[inds]
    fr = Int[1,2,3,4,5]
    try
        model = wangtest(assembly, fr; basis_translation = [1 0; 0 1; 0 0], basis_rotation = zeros(3,2))
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nBlocks 1,7 rausnehmen\n")
    inds = [1,4,5,6,7]
    assembly = ass[inds]
    fr = Int[2,3,4]
    try
        model = wangtest(assembly, fr; basis_translation = [1 0; 0 1; 0 0], basis_rotation = zeros(3,2))
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nBlocks 4 rausnehmen\n")
    inds = [4,5,6]
    assembly = ass[inds]
    fr = Int[2,3]
    try
        model = wangtest(assembly, fr; basis_translation = [1 0; 0 1; 0 0], basis_rotation = zeros(3,2))
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end

    write(io, "\nBlocks 5 rausnehmen\n")
    inds = [5,6]
    assembly = ass[inds]
    fr = Int[1]
    try
        model = wangtest(assembly, fr; basis_translation = [1 0; 0 1; 0 0], basis_rotation = zeros(3,2))
        # write(io, model)
    catch e
        showerror(io, e)
        write(io, "not interlocking. (as a result of: Error.)")
    end
    write(io, "\n------------------------------------------------\n")
end
    

close(io)
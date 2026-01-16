module Tools

using Interpolations
using NonlinearSolve
using ..DataStructures
using ..IO: init

export GridMapping
export shiftgrid!, transform_grid, calculate_grid
export search_inifile
export to_single_precision


module GridMapping
    using ..Interpolations, ForwardDiff, ..DataStructures

    export mapping, dmapping, ddmapping, dddmapping
    export spacing, stretching, dstretching, dstretching

    @inline profile(s, h1h0, δ1, st1) = (h1h0 - 1)*δ1*log(1 + exp((s - st1)/δ1))
    @inline dprofile(s, h1h0, δ1, st1) = (h1h0 - 1)/2*(1 + tanh((s - st1)/(2*δ1)))
    @inline ddprofile(s, h1h0, δ1, st1) = begin
       (h1h0 - 1)/(4*δ1)*(sech((st1 - s)/(2*δ1)))^2 
    end
    @inline dddprofile(s, h1h0, δ1, st1) = begin
        a = (st1 - s)/(2*δ1)
        (h1h0 - 1)/(4*δ1^2)*tanh(a)*(sech(a))^2
    end

    @inline mapping(s, h1h0, δ1, st1, h2h0, δ2, st2) = begin
        C = - profile(0, h1h0, δ1, st1) - profile(0, h2h0, δ2, st2)
        s + profile(s, h1h0, δ1, st1) + profile(s, h2h0, δ2, st2) + C
    end
    @inline dmapping(s, h1h0, δ1, st1, h2h0, δ2, st2) = begin
        1 + dprofile(s, h1h0, δ1, st1) + dprofile(s, h2h0, δ2, st2)
    end
    @inline ddmapping(s, h1h0, δ1, st1, h2h0, δ2, st2) = begin
        ddprofile(s, h1h0, δ1, st1) + ddprofile(s, h2h0, δ2, st2)
    end
    @inline dddmapping(s, h1h0, δ1, st1, h2h0, δ2, st2) = begin
        dddprofile(s, h1h0, δ1, st1) + dddprofile(s, h2h0, δ2, st2)
    end
    
    # NOTE Δz is also available by simple Δz=diff(z), giving the data point margins
    # But, then length(Δz) = length(z)-1
    @inline spacing(nodes, z) = begin
        itp = extrapolate(
            interpolate((nodes,), z, Gridded(Linear())), 
            Line()
        )
        return ForwardDiff.derivative.(Ref(itp), nodes)
    end
    @inline spacing(grid::Grid) = spacing(1:grid.nz, grid.z)
    @inline spacing(z::Vector{<:AbstractFloat}) = spacing(1:length(z), z)

    @inline stretching(z) = begin
        Δz = diff(z)
        (Δz[2:end] .- Δz[1:end-1])./Δz[1:end-1].*100
        # With the below variant one gets stretching in the same dimension as z
        # itp = extrapolate(
        #     interpolate((collect(1:lenth(z)),), spacing(z), Gridded(Linear())), 
        #     Line()
        # )
        # return ForwardDiff.derivative.(Ref(itp), nodes) ./ spacing(z)
    end
    @inline stretching(grid::Grid) = stretching(grid.z)
    @inline dstretching(z) = begin
        itp = extrapolate(
            interpolate((collect(1:length(z)-2),), stretching(z), Gridded(Linear())), 
            Line()
        )
        ForwardDiff.derivative.(Ref(itp), z)
    end
end


using .GridMapping


let
    """
        transform_grid(data, newgrid)
    Transform the grid of _data_ in _grid_. _shift_ corresponds to the axes given 
    in _shiftaxis_.
    """
    global function transform_grid(
            data::ScalarData,
            grid::Grid;
            shift::Vector{<:AbstractFloat}=[0.0],
            shiftaxis::Vector{Symbol}=[:z]
        )
        return _transform_grid(data, grid, shift, shiftaxis)
    end


    global function transform_grid(
            data::ScalarData,
            gridfile::String;
            shift::Vector{<:AbstractFloat}=[0.0],
            shiftaxis::Vector{Symbol}=[:z]
        )
        return _transform_grid(data, loadgrid(gridfile), shift, shiftaxis)
    end


    global function transform_grid(
            datafile::String,
            gridfile::String;
            shift::Vector{<:AbstractFloat}=[0.0],
            shiftaxis::Vector{Symbol}=[:z]
        )
        return _transform_grid(load(datafile), loadgrid(gridfile), shift, shiftaxis)
    end


    @inline function _transform_grid(
            data::ScalarData, grid::Grid, shift::Vector{<:AbstractFloat}, shiftaxis::Vector{Symbol}
        )::ScalarData
        # Container for transformed data
        newdata = init(grid)
        newdata.name = data.name
        newdata.time = data.time
        # Shift original grid
        for (i, axis) ∈ enumerate(shiftaxis)
            shiftgrid!(data, shift[i], axis=axis)
        end
        # Use interpolation of original data with higher resolution ot fill the container with the lower reolution grid
        if data.grid.ny==1
            itp = interpolate(
                (data.grid.x, data.grid.z),                 # Nodes of the grid
                data.field[:,1,:],                          # Field to be interpolated
                Gridded(Linear())                           # Interpolation type
            )
        else
            itp = interpolate(
                (data.grid.x, data.rgid.y, data.grid.z),    # Nodes of the grid
                data.field[:,:,:],                          # Field to be interpolated
                Gridded(Linear())                           # Interpolation type
            )
        end
        for k ∈ eachindex(newdata.grid.z)
            for j ∈ eachindex(newdata.grid.y)
                for i ∈ eachindex(newdata.grid.x)
                    newdata.field[i,j,k] = itp(newdata.grid.x[i], newdata.grid.z[k])
                end
            end
        end
        return newdata
    end
end


"""
    shiftgrid!(data, shift, axis)
Shift the grid of _data_ by _shift_ along _axis_, i.e. each point of _axis_ is 
added by _shift_.
"""
function shiftgrid!(data::ScalarData, shift::AbstractFloat; axis::Symbol=:z)
    if axis==:z
        data.grid.z .= data.grid.z .+ shift
    elseif axis==:y
        data.grid.y .= data.grid.y .+ shift
    elseif axis==:x
        data.grid.x .= data.grid.x .+ shift
    end
end


"""
    calculate_grid(nx, ny, nz, lx, ly, lz, maxstretching, st1) -> Grid
Returns _Grid_ with nx*ny*nz grid points, axis lengths lx, ly, and lz and 
non-uniform z-axis. The stretching along z has its peak wiht _maxstretching_. 
_st1_ represents the position of the mapping's inflection point and has a 
default value of 0.0.
"""
function calculate_grid(
        nx::Int, ny::Int, nz::Int,
        lx::AbstractFloat, ly::AbstractFloat, lz::AbstractFloat,
        maxstretching::AbstractFloat;
        st1::AbstractFloat=0.0,
        bufferlength::AbstractFloat=4*2π/abs(sin(-π/4))
    )::Grid
    """
        h0 -> intial uniform grid step
        h1 -> new maximal grid step after stretching (top of tanh function value)
        δ1 -> stretching length (tanh width)
        s -> transition to stretching occuring at s=st1

        Notes on the domain choice:
        ⋅ Because of parallelization, nz has to be a multiple of 128.
        ⋅ The grid mapping parameters are computed with Newton-Rhapson for a given 
            domain length lz and a given peak for the stretching in %
    """

    println("Calculating non-uniform grid with ")
    println("   lz=$lz, maxstretching=$maxstretching, st1=$st1, nz=$nz")
        
    # Grid points
    Nx = nx + 1; Ny = ny + 1
    Δx = lx/Nx; Δy = ly/Ny
    lz0 = (nz-1)*Δx
    st2 = lz0 - st1
    s = collect(1:nz)
    z0 = range(0.0, lz0, nz)

    # Solve the non-linear system defined by f for δ1 and h1h0
    f(u, p) = begin
        z() = mapping.(z0, u[1], u[2], st1, u[1], -u[2], st2)
        str() = stretching(z())
        [
            # Eq. for h1h0 with p[1]=lz
            mapping(lz0, u[1], u[2], st1, u[1], -u[2], st2) .- p[1],
            # Eq. for δ1 with p[2]=max.stretching
            maximum(str()) .- p[2]
        ]
    end
    problem = NonlinearProblem(
        f, 
        [100.0, -0.4], 
        [lz, maxstretching],
    )
    solution = solve(problem, SimpleNewtonRaphson())
    δ1 = solution.u[2]; δ2 = -δ1
    h1h0 = solution.u[1]; h2h0 = h1h0
    vals_1 = [st1, h1h0, δ1, st2, h2h0, δ2]
    z = mapping.(z0, h1h0, δ1, st1, h2h0, δ2, st2)

    # Print buffer info
    points = findmin(abs.(z .- bufferlength))[2]
    println("   Buffer length: $(z[points])")
    println("   Points=$points")

    printstyled("
        [BufferZone]",
        bold = false
    )
    println("
        Type=relaxation
        LoadBuffer=no
        PointsUKmax=$(points)
        PointsSKmax=$(points)
        PointsUKmin=$(points)
        PointsSKmin=$(points)
        ParametersUKmax=0.25, 3.0
        ParametersSKmax=0.25, 3.0
        ParametersUKmin=0.25, 3.0
        ParametersSKmin=0.25, 3.0
    ")

    # Print vertical grid info
    printstyled("
        [IniGridOz]", 
        bold = false
    )
    println("
        periodic=no
        segments=1

        points_1=$nz
        scales_1=$(lz0)
        opts_1=Tanh
        vals_1=$(vals_1[1]), $(vals_1[2]), $(vals_1[3]), $(vals_1[4]), $(vals_1[5]), $(vals_1[6])
    ")

    return Grid(
        nx=nx, ny=ny, nz=nz,
        lx=lx, ly=ly, lz=lz,
        x = collect(range(0.0, lx-Δx, nx)),
        y = collect(range(0.0, ly-Δy, ny)),
        z = mapping.(z0, h1h0, δ1, st1, h2h0, δ2, st2)
    )
end


function calculate_grid(
        nx::Int, ny::Int, nz::Int,
        lx::Real, ly::Real, lz::Real;
        lower_maxstretching::Real=-2.0,
        upper_maxstretching::Real=2.0,
        st1::Real=0.0,
        Δst2::Real=st1,
        lower_bufferlength::Real=4*2π/abs(sin(-π/4)),
        upper_bufferlength::Real=4*2π/abs(sin(-π/4))
    )::Grid
    """
        h0 -> intial uniform grid step
        h1 -> new maximal grid step after stretching (top of tanh function value)
        δ1 -> stretching length (tanh width)
        s -> transition to stretching occuring at s=st1

        Notes on the domain choice:
        ⋅ Because of parallelization, nz has to be a multiple of 128.
        ⋅ The grid mapping parameters are computed with Newton-Rhapson for 
            above specified constraints
    """

    println("Calculating non-uniform grid ... ")
        
    # Grid points
    Nx = nx + 1; Ny = ny + 1
    Δx = lx/Nx; Δy = ly/Ny
    lz0 = (nz-1)*Δx
    # s = collect(1:nz)
    z0 = range(0.0, lz0, nz)
    println("Δx = $Δx")
    println("Δz = $(z0[2] - z0[1])")

    if (z0[2] - z0[1])!==Δx
        error("No proper grid init")
    end

    println(lz0)

    st2 = lz0 - Δst2

    # Solve the non-linear system defined by f for δ1 and h1h0
    f(u, p) = begin
        z() = mapping.(z0, u[1], u[2], st1, u[1], u[3], st2)
        str() = stretching(z())
        [
            # Constraint for final domain length with p[1]=lz
            mapping(lz0, u[1], u[2], st1, u[1], u[3], st2) .- p[1],
            # Constraint for lower maxstretching with p[2]=lower_maxstretching
            minimum(str()) .- p[2],
            # Constraint for upper maxstretching with p[3]=upper_maxstretching
            maximum(str()) .- p[3],
        ]
    end
    problem = NonlinearProblem(
        f, 
        [200.0, -0.7, 0.7], # [h1h0, δ1, δ2] startvalues
        [lz, lower_maxstretching, upper_maxstretching], # p
    )
    solution = solve(problem, NewtonRaphson())
    h1h0 = solution.u[1]
    δ1 = solution.u[2]
    h2h0 = solution.u[1]
    δ2 = solution.u[3]
    vals_1 = [st1, h1h0, δ1, st2, h2h0, δ2]
    z = mapping.(z0, h1h0, δ1, st1, h2h0, δ2, st2)

    # Print buffer info
    lpoints = findmin(abs.(z .- lower_bufferlength))[2]
    upoints = length(z) - findmin(abs.(z .- upper_bufferlength))[2]
    println("   Lower buffer length: $(z[lpoints])")
    println("   Points=$lpoints")
    println("   Upper buffer length: $(z[end-upoints])")
    println("   Points=$upoints")

    printstyled("
        [BufferZone]",
        bold = false
    )
    println("
        Type=relaxation
        LoadBuffer=no
        PointsUKmax=$(upoints)
        PointsSKmax=$(upoints)
        PointsUKmin=$(lpoints)
        PointsSKmin=$(lpoints)
        ParametersUKmax=0.25, 3.0
        ParametersSKmax=0.25, 3.0
        ParametersUKmin=0.25, 3.0
        ParametersSKmin=0.25, 3.0
    ")

    # Print vertical grid info
    printstyled("
        [IniGridOz]", 
        bold = false
    )
    println("
        periodic=no
        segments=1

        points_1=$nz
        scales_1=$(lz0)
        opts_1=Tanh
        vals_1=$(vals_1[1]), $(vals_1[2]), $(vals_1[3]), $(vals_1[4]), $(vals_1[5]), $(vals_1[6])
    ")

    return Grid(
        nx=nx, ny=ny, nz=nz,
        lx=lx, ly=ly, lz=lz,
        x = collect(range(0.0, lx-Δx, nx)),
        y = collect(range(0.0, ly-Δy, ny)),
        z = mapping.(z0, h1h0, δ1, st1, h2h0, δ2, st2)
    )
end


"""
    search_inifile(inifile, block, key)
Search the tlab.ini-file for the value of _key_. To avoid ambigious keys, the 
_block_has to be provided, too. Blocks appar in tlab.ini in the pattern 
_[block]_.
"""
function search_inifile(file::String, block::String, key::String)::String
    f = open(file, "r")
    res = ""
    corr_block = false
    while ! eof(f)
        s = readline(f)
        # Correct block?
        if startswith(s, "[")
            if occursin(block, s)
                corr_block = true
            else
                corr_block = false
            end
        end
        # Check for key if in correct block
        if corr_block && ! startswith(s, "[")
            val = split(s, "=")
            if length(val)==2
                if key==val[1]
                    res = val[2]
                end
            end
        end
    end
    close(f)
    return res
end


"""
    to_single_precision(file1, file2)
Reads _file1_, converts the field data from Float64 to Float32 and stores 
to file2. Purpose is that field files can be used with the 
FileDataType=Single options in tlab.ini, if the data is originally saved 
in double precission.
"""
function to_single_precision(infile::String, outfile::Union{String, SubString})
    # Reading the infile
    istream = open(infile, "r")
    headersize = read(istream, Int32) # 4 bytes
    steps = Vector{Int32}(undef, 4) # step[:]=[nx,ny,nz,nt] corresponds to 4*4=16 bytes
    for i ∈ eachindex(steps)
        steps[i] = read(istream, Int32)
    end
    params = Vector{Float64}(undef, (headersize - 5*sizeof(headersize))÷sizeof(Float64))
    for i ∈ eachindex(params)
        params[i] = read(istream, eltype(params))
    end
    field = Vector{Float64}(undef, steps[1]*steps[2]*steps[3])
    read!(istream, field)
    close(istream)
    # Writing in single precission to outfile
    open(outfile, "w") do ostream
        write(ostream, headersize)
        write(ostream, steps)
        write(ostream, params)
        write(ostream, convert(Vector{Float32}, field))
    end
end


"""
    to_single_precision!(field, infile, outfile)
Muating variant of to_single_precision. Does not allocate _field_ but 
writes into the given array. Might be helpful for processing multiple files ...
"""
function to_single_precision!(
        field::Vector{<:AbstractFloat}, 
        infile::String, 
        outfile::Union{String, SubString}
    )
    # Reading the infile
    istream = open(infile, "r")
    headersize = read(istream, Int32) # 4 bytes
    steps = Vector{Int32}(undef, 4) # step[:]=[nx,ny,nz,nt] corresponds to 4*4=16 bytes
    for i ∈ eachindex(steps)
        steps[i] = read(istream, Int32)
    end
    params = Vector{Float64}(undef, (headersize - 5*sizeof(headersize))÷sizeof(Float64))
    for i ∈ eachindex(params)
        params[i] = read(istream, eltype(params))
    end
    read!(istream, field)
    close(istream)
    # Writing in single precission to outfile
    open(outfile, "w") do ostream
        write(ostream, headersize)
        write(ostream, steps)
        write(ostream, params)
        write(ostream, convert(Vector{Float32}, field))
    end
end


function to_single_precision(
        infiles::Vector{String}, 
        outfile::Union{Vector{String}, Vector{SubString}}
    )
    """
        This function preallocates buffer arrays and therefore assumes that 
        only the data changes, thus all files have the same grid.
    """
    # TODO not working ...
    # Preallocation with first infile
    istream = open(infile, "r")
    headersize = read(istream, Int32) # 4 bytes
    steps = Vector{Int32}(undef, 4) # step[:]=[nx,ny,nz,nt] corresponds to 4*4=16 bytes
    params = Vector{Float64}(undef, (headersize - 5*sizeof(headersize))÷sizeof(Float64))
    field = Vector{Float64}(undef, steps[1]*steps[2]*steps[3])
    close(istream)
    
    # Do the precision transformation
    for infile ∈ infiles
        open(infile, "r") do istream
            for i ∈ eachindex(steps)
                steps[i] = read(istream, Int32)
            end
            for i ∈ eachindex(params)
                params[i] = read(istream, eltype(params))
            end
            read!(istream, field)
        end        
        open(outfile, "w") do ostream
            write(ostream, headersize)
            write(ostream, steps)
            write(ostream, params)
            write(ostream, convert(Vector{Float32}, field))
        end
    end
end


end
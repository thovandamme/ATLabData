module Tools

using Interpolations
using SimpleNonlinearSolve
using ..DataStructures
using ..IO: init

export GridMapping
export shiftgrid!, transform_grid, calculate_grid
export search_inifile


module GridMapping
    using Interpolations, ForwardDiff, ..DataStructures

    export mapping, spacing, stretching

    @inline profile(s, h1h0, δ1, st1) = (h1h0 - 1)*δ1*log(1 + exp((s - st1)/δ1))
    @inline mapping(s, h1h0, δ1, st1, h2h0, δ2, st2) = begin
        C = - profile(0, h1h0, δ1, st1) - profile(0, h2h0, δ2, st2)
        s + profile(s, h1h0, δ1, st1) + profile(s, h2h0, δ2, st2) + C
    end
    @inline spacing(nodes, z) = begin
        itp = linear_interpolation(nodes, z, extrapolation_bc = Line())
        return ForwardDiff.derivative.(Ref(itp), nodes)
    end
    @inline spacing(grid::Grid) = begin
        return spacing(range(1:grid.nz), grid.z)
    end
    @inline spacing(z) = spacing(range(1:length(z)), z)
    @inline stretching(nodes, z) = begin
        itp = linear_interpolation(nodes, spacing(nodes, z), extrapolation_bc = Line())
        return ForwardDiff.derivative.(Ref(itp), nodes) ./ spacing(nodes, z) .* 100
    end
    @inline stretching(grid::Grid) = begin
        return stretching(reange(1:grid.nz), grid.z)
    end
    @inline stretching(z) = stretching(range(a:length(z)), z)
end


using .GridMapping


let
    """
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
Shift the grid of _data_ by _shift along _axis_, i.e. each point of _axis_ is 
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
    lz0 = nz*lx/nx
    st2 = lz0 - st1
    s = collect(1:nz)
    z0 = range(0.0, lz0, nz)

    # Solve the non-linear system defined by f for δ1 and h1h0
    f(u, p) = [
        # Eq. for h1h0 with p[1]=lz
        mapping(lz0, u[1], u[2], st1, u[1], -u[2], st2) .- p[1],
        # Eq. for δ1 with p[2]=max.stretching
        maximum(stretching(s, mapping.(z0, u[1], u[2], st1, u[1], -u[2], st2))) .- p[2],
    ]
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
        scales_1=$(nz*lx/nx)
        opts_1=Tanh
        vals_1=$(vals_1[1]), $(vals_1[2]), $(vals_1[3]), $(vals_1[4]), $(vals_1[5]), $(vals_1[6])
    ")

    # # Print info about the spacing error
    # δz = 2*Δz
    # imin = findmin(abs.(z .- (zs - δz)))[2]
    # imax = findmin(abs.(z .- (zs + δz)))[2]
    # err_imin = spacing[imin] - lx/nx
    # err_imax = spacing[imax] - lx/nx
    # printstyled("\n Spacing deviation for nz=$nz, st1=$st1, stretching=$stretch: \n", bold=true)
    # if err_imin > 1.e-4 || err_imax > 1.e-4
    #     printstyled("  $err_imin  at z=$(z[imin]) \n", color=:light_red)
    #     printstyled("  $err_imax  at z=$(z[imax]) \n", color=:light_red)
    # else
    #     printstyled("  $err_imin  at z=$(z[imin]) \n", color=:cyan)
    #     printstyled("  $err_imax  at z=$(z[imax]) \n", color=:cyan)
    # end

    return Grid(
        nx=nx, ny=ny, nz=nz,
        lx=lx, ly=ly, lz=lz,
        x = collect(range(0.0, lx, nx)),
        y = collect(range(0.0, ly, ny)),
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


end
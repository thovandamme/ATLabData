module IO

using NetCDF
using ..DataStructures

export load, load!, loadgrid
export init


# ------------------------------------------------------------------------------
#                           API (exported functions)
# ------------------------------------------------------------------------------
"""
    loadgrid(gridfile) -> Grid
Loads the grid data from the file _gridfile_ into the composite type _grid_.
"""
loadgrid(gridfile::String)::Grid = Grid_from_gridfile(gridfile)


"""
    load(file) -> ScalarData
Load the data contained in the path _file_ into the type _Data_.

    load(dir, field, time) -> ScalarData
Search in _dir_ for the file containing _field_ at _time_ and load into the type 
_Data_.

    load(xfile, yfile, zfile) -> VectorData
Load the data contained in the paths _xfile_, _yfile_ and _zfile_ into the 
type _VectorData_.

    load(file, var) -> AveragesData
Load the data contained in the path _file_ into the type _AveragesData_.
_file_ has to be NetCDF file containing the averages from _average.x_.
"""
load(file::String)::ScalarData = ScalarData_from_file(file)
load(dir::String, field::String, time::Real; component::String=".0")::ScalarData = load(file_for_time(dir, time, field, component))

load(xfile::String, yfile::String, zfile::String)::VectorData = VectorData_from_files(xfile, yfile, zfile)

load(file::String, var::String)::AveragesData = AveragesData_from_NetCDF(file, var)


"""
    load!(data, file)
Version of _load_ for preallocated data container. This function does not 
    update the grid attribute!
"""
load!(data::ScalarData, file::String) = ScalarData_from_file!(data, file)
load!(
    data::VectorData, 
    xfile::String, yfile::String, zfile::String
) = VectorData_from_file!(data, xfile, yfile, zfile)


"""
    init(grid; dims=1) -> AbstractData
Initialize an empty data container. dims has to be an integer in of value 
1 or 3. 1 correpsonds to _ScalarData_ while 3 returns _VectorData_. 
Default for _dims_ is 1. _grid_ has to be of type _Grid_.
"""
function init(
        grid::Grid{T,I};
        dims::Int=1
    )::AbstractData where {T<:AbstractFloat, I<:Signed}
    if dims==1
        return ScalarData(
            name = "empty container",
            grid = grid,
            time = convert(T, 0.0),
            field = Array{T, 3}(undef, grid.nx, grid.ny, grid.nz)
        )
    elseif dims==3
        return VectorData(
            name = "empty container",
            grid = grid,
            time = convert(T, 0.0),
            xfield = Array{T, 3}(undef, grid.nx, grid.ny, grid.nz),
            yfield = Array{T, 3}(undef, grid.nx, grid.ny, grid.nz),
            zfield = Array{T, 3}(undef, grid.nx, grid.ny, grid.nz),
        )
    end
end


function init(gridfile::String; dims::Int=1, T::Type=Float64)::AbstractData
    return init(convert(T, loadgrid(gridfile)), dims=dims)
end


# ------------------------------------------------------------------------------
#                                   Grid
# ------------------------------------------------------------------------------
""" 
    Load the file _grid_ into composite type _Grid_
"""
function Grid_from_gridfile(gridfile::String)::Grid
    "f90_record_markers from the amount of write commands in storing routine "
    io = open(gridfile, "r")
    marker1 = read(io, Int32)
    nx = read(io, Int32)
    ny = read(io, Int32)
    nz = read(io, Int32)
    marker2 = read(io, Int32)
    marker3 = read(io, Int32)
    scalex = read(io, Float64)
    scaley = read(io, Float64)
    scalez = read(io, Float64)
    marker4 = read(io, Int32)
    marker5 = read(io, Int32)
    x = Vector{Float64}(undef, nx)
    for i ∈ 1:nx
        x[i] = read(io, Float64)
    end
    marker7 = read(io, Int32)
    marker8 = read(io, Int32)
    y = Vector{Float64}(undef, ny)
    for i ∈ 1:ny
        y[i] = read(io, Float64)
    end
    marker9 = read(io, Int32)
    marker10 = read(io, Int32)
    z = Vector{Float64}(undef, nz)
    for i ∈ 1:nz
        z[i] = read(io, Float64)
    end
    return Grid(nx, ny, nz, scalex, scaley, scalez, x, y, z)
end


function Grid_from_file(dir::String)::Grid
    gridfile = joinpath(dir, "grid")
    return Grid_from_gridfile(gridfile)
end


# ------------------------------------------------------------------------------
#                               ScalarData
# ------------------------------------------------------------------------------
"""
    Loading data that is stored in a single file.
"""
function ScalarData_from_file(fieldfile::String)::ScalarData
    verbose("ScalarData", fieldfile)
    filename = split(fieldfile, "/")[end]
    if startswith(filename, "flow.") || startswith(filename, "scal.")
        return ScalarData_from_raw(fieldfile)
    else
        return ScalarData_from_visuals(fieldfile)
    end
end


function ScalarData_from_raw(fieldfile::String)::ScalarData{Float64, Int32}
    grid = Grid_from_file(dirname(fieldfile))
    buffer, t = Array_from_rawfile(grid, fieldfile)
    return ScalarData(
        name = splitpath(fieldfile)[end],
        grid = grid,
        time = t,
        field = buffer
    )
end


function ScalarData_from_visuals(fieldfile::String)::ScalarData{Float32, Int32}
    grid = convert(Float32, Grid_from_file(dirname(fieldfile)))
    return ScalarData(
        name = splitpath(fieldfile)[end],
        grid = grid,
        time = time_from_file(fieldfile),
        field = Array_from_file(grid, fieldfile)
    )
end


function ScalarData_from_file!(data::ScalarData, fieldfile::String)
    filename = split(fieldfile, "/")[end]
    if eltype(data)[1]==Float64 && !(startswith(filename, "flow.") || startswith(filename, "scal."))
        error(
            "_data_ is initialized with type $(eltype(data)[1]) but $filename 
            contains floats of type Float32"
        )
    elseif eltype(data)[1]==Float32 && (startswith(filename, "flow.") || startswith(filename, "scal."))
        error(
            "_data_ is initialized with type $(eltype(data)[1]) but $filename 
            contains floats of type Float64"
        )
    end
    verbose("ScalarData", fieldfile)
    data.name = filename
    if startswith(filename, "flow.") || startswith(filename, "scal.")
        data.field, data.time = Array_from_rawfile(data.grid, fieldfile)
    else
        data.time = time_from_file(fieldfile)
        data.field .= Array_from_file(data.grid, fieldfile)
    end
    return nothing
end


# ------------------------------------------------------------------------------
#                               VectorData
# ------------------------------------------------------------------------------
"""
    Loading data from three files. No check for physical consistency.
"""
function VectorData_from_files(
        xfieldfile::String,
        yfieldfile::String,
        zfieldfile::String
    )::VectorData
    verbose("VectorData", xfieldfile, yfieldfile, zfieldfile)
    filename = split(xfieldfile, "/")[end]
    if startswith(filename, "flow.") || startswith(filename, "scal.")
        return VectorData_from_raw(xfieldfile, yfieldfile, zfieldfile)
    else
        return VectorData_from_visuals(xfieldfile, yfieldfile, zfieldfile)
    end
end


function VectorData_from_raw(
        xfieldfile::String,
        yfieldfile::String,
        zfieldfile::String
    )::VectorData{Float64, Int32}
    grid = Grid_from_file(dirname(xfieldfile))
    buffer, t = Array_from_rawfile(grid, xfieldfile)
    return VectorData(
        name = string(splitpath(xfieldfile)[end][1:end-2]),
        grid = grid,
        time = t,
        xfield = buffer,
        yfield = Array_from_rawfile(grid, yfieldfile)[1],
        zfield = Array_from_rawfile(grid, zfieldfile)[1]
    )
end


function VectorData_from_visuals(
        xfieldfile::String,
        yfieldfile::String,
        zfieldfile::String
    )
    grid = convert(Float32, Grid_from_file(dirname(xfieldfile)))
    return VectorData(
        name = string(splitpath(xfieldfile)[end][1:end-2]),
        grid = grid,
        time = time_from_file(xfieldfile),
        xfield = Array_from_file(grid, xfieldfile),
        yfield = Array_from_file(grid, yfieldfile),
        zfield = Array_from_file(grid, zfieldfile)
    )
end


function VectorData_from_files!(
        data::VectorData, 
        xfieldfile::String,
        yfieldfile::String,
        zfieldfile::String
    )
    filename = split(xfieldfile, "/")[end]
    if eltype(data)[1]==Float64 && !(startswith(filename, "flow.") || startswith(filename, "scal."))
        error(
            "_data_ is initialized with type $(eltype(data)[1]) but $filename 
            contains floats of type Float32"
        )
    elseif eltype(data)[1]==Float32 && (startswith(filename, "flow.") || startswith(filename, "scal."))
        error(
            "_data_ is initialized with type $(eltype(data)[1]) but $filename 
            contains floats of type Float64"
        )
    end
    verbose("VectorData", xfieldfile, yfieldfile, zfieldfile)
    data.name = string(splitpath(xfieldfile)[end][1:end-2])
    if startswith(filename, "flow.") || startswith(filename, "scal.")
        data.xfield, data.time = Array_from_rawfile(data.grid, xfieldfile)
        data.yfield .= Array_from_rawfile(data.grid, yfieldfile)[1]
        data.zfield .= Array_from_rawfile(data.grid, zfieldfile)[1]
    else
        data.time = time_from_file(fieldfile)
        data.xfield .= Array_from_file(data.grid, xfieldfile)
        data.yfield .= Array_from_file(data.grid, yfieldfile)
        data.zfield .= Array_from_file(data.grid, zfieldfile)
    end
    return nothing
end


# ------------------------------------------------------------------------------
#                         Array from binary
# ------------------------------------------------------------------------------
"""
    Load array from a binary file according to the information in _grid_.
"""
function Array_from_file(
        grid::Grid{T,I}, fieldfile::String
    )::Array{T,3} where {T<:AbstractFloat, I<:Signed}
    buffer = Vector{T}(undef, grid.nx*grid.ny*grid.nz)
    read!(fieldfile, buffer)
    return reshape(buffer, (grid.nx, grid.ny, grid.nz))
end
 

function Array_from_rawfile(
        grid::Grid{T,I}, fieldfile::String
    )::Tuple{Array{T, 3}, Float64} where {T<:AbstractFloat, I<:Signed}
    io = open(fieldfile, "r")
    """ 
        Header contains: 
        headersize/offset, nx, ny, nz, nt, time, visc, froude, schmidt
        unformatted with no record markers
    """
    headersize = read(io, Int32) # First entry is the header length in bytes
    seek(io, 5*sizeof(headersize))
    time = read(io, Float64)
    seek(io, headersize) # Jump to first last entry belonging to the header
    buffer = Vector{Float64}(undef, grid.nx*grid.ny*grid.nz)
    read!(io, buffer)
    return reshape(buffer, (grid.nx, grid.ny, grid.nz)), time
end


# ------------------------------------------------------------------------------
#                            AveragesData
# ------------------------------------------------------------------------------
function AveragesData_from_NetCDF(file::String, var::String)::AveragesData
    return AveragesData(
        name = var,
        time = ncread(file, "t"),
        z = ncread(file, "z"),
        field = ncread(file, var)
    )
end


# ------------------------------------------------------------------------------
#                             Helping functions
# ------------------------------------------------------------------------------
function timestep_from_filename(filename::String)::Int
    namestring = split(filename, "/")[end]
    namestring = split(namestring, ".")[1]
    stepstring = namestring[end-5:end]
    return parse(Int, stepstring)
end


function time_from_timestep(timestep::Int, dir::String)::AbstractFloat
    avgfile = "avg"*string(timestep)*".nc"
    avgfile = joinpath(dir, avgfile)
    return ncread(avgfile, "t")[1]
end


function time_from_file(file::String)::Float32
    tstep = timestep_from_filename(file)
    t = time_from_timestep(tstep, dirname(file))
    return t
end


function file_for_time(
        dir::String, 
        time::Real,
        field::String,
        component::String=".0"
    )::String
    Δt = time
    f = ""
    maxΔt_relative = 0.1
    maxΔt = maxΔt_relative*time
    println("Searching for t = ", time)
    for filename in readdir(dir)
        correct_field = false
        if component ∈ (".1", ".2", ".3", ".4", ".5", ".6")
            if startswith(filename, field) && endswith(filename, component)
                correct_field = true
            end
        elseif component == ".0"
            if startswith(filename, field)
                correct_field = true
            end
        else
            error("Loading a field with this _component_ is not implemented.")
        end
        if correct_field
            file = joinpath(dir, filename)
            t = time_from_file(file)
            printstyled("   ", filename, "    ", t, "\n", color=:cyan)
            if abs(time - t) < Δt
                f = file
                Δt = abs(time - t)
            end
        end
    end
    if Δt > maxΔt
        @warn "The time loaded exceeds the maximal accepted relative error of $(maxΔt_relative)."
    end
    return f
end


function verbose(dtype::String, files...; save::Bool=false)
    if !save
        text = "Loading " * dtype * " ..."
    else
        text = "Saving " * dtype * " ..."
    end
    println(text)
    for file in files
        printstyled("   "*file*"\n", color=:cyan)
    end
end


function done()
    printstyled("   Done \n", color=:green)
end


end
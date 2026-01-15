module Basics


import Base: +, -, *, /, ^, size, display, abs, convert, deepcopy, eltype

using ..DataStructures

export size, display, +, -, *, /, ^, abs, log, view, convert, eltype
export crop, logarithm
export component
export rms, average, mean, flucs


"""
    size(data)
Returns `sìze(data.field)`.
"""
size(data::ScalarData)::Tuple = size(data.field)
size(data::VectorData)::Tuple = (size(data.xfield), size(data.yfield), size(data.zfield))
size(data::AveragesData)::Tuple = size(data.field)


"""
    data1 - data2
Subtract the _field_ of _data2_ from the _field_ of _data1_. 
Does the same as '''subtract(data1, data2)'''
"""
-(data1::ScalarData, data2::ScalarData)::ScalarData = ScalarData(
    name = "$(data1.name)-$(data2.name)",
    field = data1.field .- data2.field,
    time = data1.time,
    grid = data1.grid
)
-(data1::VectorData, data2::VectorData)::VectorData = VectorData(
    name = "$(data1.name)-$(data2.name)",
    grid = data1.grid,
    time = data1.time,
    xfield = data1.xfield .- data2.xfield,
    yfield = data1.yfield .- data2.yfield,
    zfield = data1.zfield .- data2.zfield
)
-(data1::AveragesData, data2::AveragesData)::AveragesData = AveragesData(
    name = "$(data1.name)-$(data2.name)",
    time = data1.time,
    z = data1.grid.z,
    field = data1.field .- data2.field
)


"""
    data1 + data2
Add the _fields_ of _data1_ and _data2_. Does the same as '''add(data1, data2)'''
"""
+(data1::ScalarData, data2::ScalarData)::ScalarData = ScalarData(
    field = data1.field .+ data2.field,
    time = data1.time,
    grid = data1.grid,
    name = data1.name*"+"*data2.name
)
+(data1::VectorData, data2::VectorData)::VectorData = VectorData(
    name = data1.name * "+" * data2.name,
    grid = data1.grid,
    time = data1.time,
    xfield = data1.xfield .+ data2.xfield,
    yfield = data1.yfield .+ data2.yfield,
    zfield = data1.zfield .+ data2.zfield
)
+(data1::AveragesData, data2::AveragesData)::AveragesData = AveragesData(
    name = "$(data1.name)+$(data2.name)",
    time = data1.time,
    z = data1.grid.z,
    field = data1.field .+ data2.field
)


"""
    data * factor
Multiply the _field_ of _data_ with _factor_. Does the same as _rescale_.

    data1 * data2
Vectorized multiplication of the fields of _data1_ and _data2_
"""
*(data::ScalarData, factor::Real)::ScalarData = ScalarData(
    name = data.name,
    time = data.time,
    grid = data.grid,
    field = data.field .* convert(Float32, factor)
)
*(data::VectorData, factor::Real)::VectorData = VectorData(
    name = data.name,
    time = data.time,
    grid = data.grid,
    xfield = data.xfield * convert(Float32, factor),
    yfield = data.yfield * convert(Float32, factor),
    zfield = data.zfield * convert(Float32, factor)
)
*(data::AveragesData, factor::Real)::AveragesData = AveragesData(
    name = data.name,
    time = data.time,
    z = data.grid.z,
    field = data.field .* convert(Float32, factor)
)
*(data1::ScalarData, data2::ScalarData)::ScalarData = ScalarData(
    name = "$(data1.name)*$(data2.name)",
    time = data1.time,
    grid = data1.grid,
    field = data1.field .* data2.field
)
*(data1::AveragesData, data2::AveragesData) = AveragesData(
    name = "$(data1.name)*$(data2.name)",
    time = data1.time,
    z = data1.grid.z,
    field = data1.field .* data2.field
)
*(factor::Real, data::AbstractData)::AbstractData = data*factor


"""
    data / number
Divide the _field_ of _data_ with _number__.

    data1 / data2
Vectorized division of the fields of _data1_ and _data2_.
"""
/(data::ScalarData, factor::Real)::ScalarData = ScalarData(
    name = data.name,
    time = data.time,
    grid = data.grid,
    field = data.field ./ convert(Float32, factor)
)
/(data::VectorData, factor::Real)::VectorData = VectorData(
    name = data.name,
    time = data.time,
    grid = data.grid,
    xfield = data.xfield ./ convert(Float32, factor),
    yfield = data.yfield ./ convert(Float32, factor),
    zfield = data.zfield ./ convert(Float32, factor)
)
/(data::AveragesData, factor::Real)::AveragesData = AveragesData(
    name = data.name,
    time = data.time,
    z = data.grid.z,
    field = data.field ./ convert(Float32, factor)
)
/(data1::ScalarData, data2::ScalarData)::ScalarData = ScalarData(
    name = "$(data1.name)*$(data2.name)",
    time = data1.time,
    grid = data1.grid,
    field = data1.field ./ data2.field
)
/(data1::AveragesData, data2::AveragesData) = AveragesData(
    name = "$(data1.name)*$(data2.name)",
    time = data1.time,
    z = data1.grid,
    field = data1.field ./ data2.field
)


"""
    data^exponent
Exponentiation of the _Data_ type. Returns _Data_  with data.field.^2 while 
maintaining the metadata.
"""
^(data::ScalarData, exponent::Real)::ScalarData = ScalarData(
    name = data.name, 
    time = data.time, 
    grid = data.grid,
    field = data.field.^exponent
)
^(data::AveragesData, exponent::Real)::AveragesData = AveragesData(
    name = data.name, 
    time = data.time, 
    z = data.grid.z,
    field = data.field.^exponent
)


"""
    abs(data)
Returns absolute values of _data_ while preserving the metadata. 
For _VectorData_ the euclidian norm is return.
"""
abs(data::ScalarData) = ScalarData(
    field = abs.(data.field),
    grid = data.grid,
    name = "abs(" * data.name * ")",
    time = data.time
)
abs(data::VectorData) = norm(data)
norm(data::VectorData)::ScalarData = ScalarData(
    field = sqrt.(data.xfield.^2 .+ data.yfield.^2 .+ data.zfield.^2),
    grid = data.grid,
    name = "abs(" * data.name * ")",
    time = data.time
)


"""
    log(data)
Returns the logarithm while preserving the metadata. Where the result is Inf 
    or -Inf the value is replaced by zero.
"""
log(b::Base.Complex, data::ScalarData) = ScalarData(
    name = "log("*data.name*")",
    time = data.time,
    grid = data.grid,
    field = replace(log.(b, data.field), Inf=>0.0, -Inf=>0.0)
)
log(data::ScalarData) = ScalarData(
    name = "log("*data.name*")",
    time = data.time,
    grid = data.grid,
    field = replace(log.(data.field), Inf=>0.0, -Inf=>0.0)
)


sqrt(data::ScalarData) = ScalarData(
    name = "√($(data.name))",
    time = data.time,
    grid = data.grid,
    field = sqrt.(data.field)
)


# maximum(data::AveragesData{T,I}) where {T<:AbstractFloat} = maximum(data.field)
# maximum(data::ScalarData{T,I}) where {T<:AbstractFloat} = maximum(data.field)


# minimum(data::AveragesData{T,I}) where {T<:AbstractFloat} = minimum(data.field)
# minimum(data::ScalarData{T,I}) where {T<:AbstractFloat} = minimum(data.field)


"""
    convert(T, data)
Convert the entries in data to T.
"""
convert(T::Type{<:AbstractFloat}, grid::Grid)::Grid{T,Int32} = Grid{T,Int32}(
    grid.nx, grid.ny, grid.nz, 
    grid.scalex, grid.scaley, grid.scalez, 
    grid.x, grid.y, grid.z
)
convert(T::Type{<:AbstractFloat}, data::ScalarData)::ScalarData{T,Int32} = ScalarData{T,Int32}(
    data.name, convert(T, data.grid), data.time, data.field
)
convert(T::Type{<:AbstractFloat}, data::VectorData)::VectorData{T,Int32} = VectorData{T,Int32}(
    data.name, convert(T, data.grid), data.time, data.xfield, data.yfield, data.zfield
)


"""
    eltype(data)
"""
eltype(grid::Grid)::Tuple = (eltype(grid.x), eltype(grid.nx))
eltype(data::AbstractData)::Tuple = eltype(data.grid)


# """
#     deepcopy(data, newdata)
# Copy _data_ into _newdata_.
# """
# # deepcopy(a::ScalarData, b::ScalarData) = 



"""
    display(data)
Show the available atributes of _data_.
"""
function display(data::Grid)
    println(typeof(data), " with attributes: ")
    print("   nx: "); println(data.nx)
    print("   ny: "); println(data.ny)
    print("   nz: "); println(data.nz)
    print("   scalex: "); println(data.scalex)
    print("   scaley: "); println(data.scaley)
    print("   scalez: "); println(data.scalez)
    print("   x: "); println(typeof(data.x))
    print("   y: "); println(typeof(data.y))
    print("   z: "); println(typeof(data.z))
end
function display(data::ScalarData)
    println(typeof(data), " with attributes: ")
    print("   name: "); println(data.name)
    print("   grid: "); println(typeof(data.grid))
    print("   field: "); println(typeof(data.field))
    print("   time: "); println(data.time)
    println("   Use display(data.grid) and display(data.field) to print the corresponding attributes")
end
function display(data::VectorData)
    println(typeof(data), " with attributes: ")
    print("   name: "); println(data.name)
    print("   grid: "); println(typeof(data.grid))
    print("   xfield: "); println(typeof(data.xfield))
    print("   yfield: "); println(typeof(data.yfield))
    print("   zfield: "); println(typeof(data.zfield))
    print("   time: "); println(data.time)
    println("   Use display(data.grid) and display(data.*field) to print the corresponding attributes")
end
function display(data::AveragesData)
    println(typeof(data), " with attributes: ")
    print("   name: "); println(data.name)
    print("   z: "); println(typeof(data.grid.z))
    print("   field: "); println(typeof(data.field))
    print("   time: "); println(data.time)
end


# ------------------------------------------------------------------------------
# -------------------- Additional operations -----------------------------------
# ------------------------------------------------------------------------------
"""
    component(data, field) -> ScalarData
Extract from data of type _VectorData_ the component _field_ 
while keeping the metadata.
"""
function component(data::VectorData, field::Symbol)::ScalarData 
    if field == :x
        return ScalarData(
            name = data.name * " - "*string(field),
            time = data.time,
            grid = data.grid,
            field = data.xfield
        )
    elseif field == :y
        return ScalarData(
            name = data.name * " - "*string(field),
            time = data.time,
            grid = data.grid,
            field = data.yfield
        )
    elseif field == :z
        return ScalarData(
            name = data.name * " - "*string(field),
            time = data.time,
            grid = data.grid,
            field = data.zfield
        )
    else
        error("Choose the x, y, or z component of the VectorData")
    end
end


function crop(
        data::ScalarData; 
        xmin = data.grid.x[1], xmax = data.grid.x[end], 
        ymin = data.grid.y[1], ymax = data.grid.y[end],
        zmin = data.grid.z[1], zmax = data.grid.z[end]
    )::ScalarData
    println("Croping ...")
    printstyled("   "*data.name, "\n", color=:cyan)
    kmin::Int=1; kmax::Int=data.grid.nz
    jmin::Int=1; jmax::Int=data.grid.ny
    imin::Int=1; imax::Int=data.grid.nx
    inrange = false; was_inrange = false
    for k in eachindex(data.field[1,1,:])
        z = data.grid.z[k]
        if zmin <= z <= zmax
            inrange = true
        else
            inrange = false
        end
        if inrange && !(was_inrange)
            kmin = k
        end
        if !(inrange) && was_inrange
            kmax = k-1
        end
        if inrange
            was_inrange = true
        else
            was_inrange = false
        end
    end
    inrange = false; was_inrange = false
    for j in eachindex(data.field[1,:,1])
        y = data.grid.y[j]
        if ymin <= y <= ymax
            inrange = true
        else
            inrange = false
        end
        if inrange && !(was_inrange)
            jmin = j
        end
        if !(inrange) && was_inrange
            jmax = j-1
        end
        if inrange
            was_inrange = true
        else
            was_inrange = false
        end
    end
    inrange = false; was_inrange = false
    for i in eachindex(data.field[:,1,1])
        x = data.grid.x[i]
        if xmin <= x <= xmax
            inrange = true
        else
            inrange = false
        end
        if inrange && !(was_inrange)
            imin = i
        end
        if !(inrange) && was_inrange
            imax = i-1
        end
        if inrange
            was_inrange = true
        else
            was_inrange = false
        end
    end
    return ScalarData{eltype(data)[1], eltype(data)[2]}(
        data.name*"-crop",
        Grid{eltype(data)[1], eltype(data)[2]}(
            imax + 1 - imin,
            jmax + 1 - jmin,
            kmax + 1 - kmin,
            xmax - xmin,
            ymax - ymin,
            zmax - zmin,
            data.grid.x[imin:imax],
            data.grid.y[jmin:jmax],
            data.grid.z[kmin:kmax]
        ), 
        data.time,
        data.field[imin:imax,jmin:jmax,kmin:kmax],
    )
end


function _crop(
        data::ScalarData;
        xmin = data.grid.x[1], xmax = data.grid.x[end], 
        ymin = data.grid.y[1], ymax = data.grid.y[end],
        zmin = data.grid.z[1], zmax = data.grid.z[end]
    )::ScalarData
    println("Croping ...")
    printstyled("   "*data.name, "\n", color=:cyan)
    imin = findmin(abs.(data.grid.x .- xmin))[2]
    imax = findmin(abs.(data.grid.x .- xmax))[2]
    jmin = findmin(abs.(data.grid.y .- ymin))[2]
    jmax = findmin(abs.(data.grid.y .- ymax))[2]
    kmin = findmin(abs.(data.grid.k .- zmin))[2]
    kmax = findmin(abs.(data.grid.k .- zmax))[2]
    return ScalarData(
        "crop($(data.name))",
        Grid{eltype(data)[1], eltype(data)[2]}(
            imax + 1 - imin,
            jmax + 1 - jmin,
            kmax + 1 - kmin,
            data.grid.x[imax+1] - data.grid.x[imin],
            data.grid.y[jmax+1] - data.grid.y[jmin],
            data.grid.z[kmax] - data.grid.z[kmin],
            data.grid.x[imin:imax],
            data.grid.y[jmin:jmax],
            data.grid.z[kmin:kmax]
        ), 
        data.time,
        data.field[imin:imax,jmin:jmax,kmin:kmax],
    )
end


end
module Physics

using ..DataStructures, ..IO, ..Basics, ..Statistics

using ForwardDiff

export vorticity, enstrophy, Ri, tke, TurbulenceScales


"""
    vorticity(data) -> VectorData
Return the curl of data, thus is a physical alternatice to curl if data 
is a velocity field.

    vorticity(dir, time) -> VectorData
Looks for the proper velocity files in dir that are nearest to _time_ and 
computes the curl.
"""
vorticity(u::VectorData)::VectorData = curl(u)
vorticity(dir::String, time::Real)::VectorData = curl(load(
    file_for_time(dir, "VelocityVector", time, ".1"),
    file_for_time(dir, "VelocityVector", time, ".2"),
    file_for_time(dir, "VelocityVector", time, ".3")
))


"""
    enstrophy(u) -> ScalarData
Calculates the enstrophy of the givem velocity field.

    enstrophy(dir, time) -> ScalarData
Looks in _dir_ for the velocity field at _time_ and calculates the appropriate 
    enstrophy.
"""
enstrophy(u::VectorData)::ScalarData = ScalarData(
    name = "enstrophy(" * u.name * ")", 
    time = u.time, 
    grid = u.grid, 
    field = norm(vorticity(u)).field.^2
)
enstrophy(dir::String, time::Real)::ScalarData = enstrophy(load(
    file_for_time(dir, "VelocityVector", time, ".1"),
    file_for_time(dir, "VelocityVector", time, ".2"),
    file_for_time(dir, "VelocityVector", time, ".3")
))


"""
Computes and returns the local Richardson number field from the given 
buoyancy and velocity fields.

    Ri(b, u) -> ScalarData
_u_ is given as _VectorData_.

    RI(b, ux, uy, uz) -> ScalarData
The single components are given as _Data_.

    Ri(dir, time) -> ScalarData
Looks in _dir_ for the buoyancy and velocity fields for _time_.
"""
Ri(b::ScalarData, u::VectorData)::ScalarData = ScalarData(
    name = "Rig("*u.name*")",
    time = b.time,
    grid = b.grid,
    field = norm(gradient(u)).field.^2 ./ norm(gradient(b))
)
Ri(b::ScalarData, ux::ScalarData, uy::ScalarData, uz::ScalarData)::ScalarData = ScalarData(
    name = "Rig("*b.name*")",
    time = b.time,
    grid = b.grid,
    field = norm(gradient(b)).field ./ (norm(gradient(ux)).field.^2 + norm(gradient(uy)).field.^2 + norm(gradient(uz)).field.^2)
)
Ri(dir::String, time::Real)::ScalarData = Ri(
    load(dir, "Buoyancy", time),
    load(dir, "VelocityVector", time, ".1"),
    load(dir, "VelocityVector", time, ".2"),
    load(dir, "VelocityVector", time, ".3"),
)


function Reynolds_stress(u::VectorData)::Matrix
    # TODO
end


function tke(data::VectorData)::ScalarData
    buffer = flucs(data)
    return ScalarData(
        "tke($(data.name))", 
        data.grid, 
        data.time,
        0.5f0 .* (buffer.xfield.^2 .+ buffer.yfield.^2 .+ buffer.zfield.^2)
    )
end


function turbulent_diffusivity(
        flux::AveragesData, 
        mean::AveragesData;
        axis::Vector{<:AbstractFloat} = mean.grid.z
    )::AveragesData
    return AveragesData(
        name = "turbDiff($(flux.name))",
        time = flux.time,
        z = axis,
        field = turbulent_diffusivity(flux.field, mean.field, axis=axis)
    )
end


function turbulent_diffusivity(
        flux::Vector{T}, mean::Vector{T}; axis::Vector{T}
    )::Vector{T} where {T<:AbstractFloat}
    return - flux ./ ForwardDiff.derivative.(mean, axis)
end


"""
    vertical_flux(w, f) -> ⟨wf⟩
Returns the vertical flux of _f_ with _w_ as the vertical velocity component.
"""
function vertical_flux(
        w::ScalarData{T,I}, f::ScalarData{T,I}
    )::ScalaData{T,I} where {T<:AbstractFloat, I<:Signed}
    return average(w*f)
end


# module TurbulenceScales
#     time(tke::AveragesData, ε::AveragesData)::Real = maximum(tke/ε)
#     length(tke::AveragesData, ε::AveragesData)::Real = maximum(tke^(3/2)/ε)
#     velocity(tke::AveragesData)::Real = maximum(sqrt(2*tke))
# end


end
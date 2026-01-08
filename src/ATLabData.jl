"""
Module to load and process data from ATLab in Julia.
"""
module ATLabData

include("DataStructures.jl")
using .DataStructures
    export Grid, ScalarData, VectorData, AveragesData

include("IO.jl")
using .IO
    export load, load!, loadgrid
    export init
    # To be removed from exports list
    export file_for_time, Grid_from_file, VAR

include("Basics.jl")
using .Basics
    export size, display, +, -, *, ^, abs, log, convert, eltype
    export crop, norm, logarithm
    export component

include("Analysis.jl")
using .Analysis
    export gradient, curl

include("Statistics.jl")
using .Statistics
    export average, rms, mean, mean!
    export flucs, flucs!, wave, wave!, turbulence, turbulence!

include("Tools.jl")
using .Tools
    export shiftgrid!, transform_grid, calculate_grid, GridMapping
    export search_inifile
    export to_single_precision

include("Physics.jl")
using .Physics
    export vorticity, enstrophy, Ri, tke, TurbulenceScales


function __init__()
    # 
end


end
using ATLabData
using Test

@testset "ATLabData.jl" begin
    
    grid = loadgrid("/home/thomas/simulations/tmp/WTD2025/dU2.0/grid")
    data = convert(Float32, init(grid))
    load!(data, "/home/thomas/simulations/tmp/WTD2025/dU2.0/VorticityVector019200.2")
end

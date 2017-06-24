using RemPiO2
using Base.Test


# write your own tests here

x = 0.32423423
out1 = RemPiO2.rem_pio2(x)
out2 = Base.Math.ieee754_rem_pio2(x)

@test out1[1] == out2[1]
@test out1[2].hi == out2[2][1]
@test out1[2].lo == out2[2][2]

using BenchmarkTools
inital_bound = 0.0
intervals = []
# Testing some important intervals and values around bad cases that need higher precision
for i = 1:7
    push!(intervals, linspace(-pi*i/4, -pi*(i-1)/4, 100000))
    push!(intervals, linspace(pi*(i-1)/4, pi*i/4, 100000))
end

push!(intervals, linspace(-2.0^20, -pi*7/4, 100000))
push!(intervals, linspace(pi*7/4, 2.0^20, 100000))
push!(intervals, linspace(2.0^20, 2.0^60, 100000))
#intervals = [intervals[end]]
time_my_rem = []
time_rem = []
println()
println("Testing speed and accuracy of rempio2")
println()
println("(1) -> New rempio2")
println("(2) -> Old rempio2")
println("(3) -> (1)/(2)")
println()
println("(1)        (2)        (3)      ")
println("--"^17)

for (i, inter) in enumerate(intervals)
    x = collect(inter)
    push!(time_my_rem, @belapsed RemPiO2.rem_pio2.(x))
    push!(time_rem, @belapsed Base.Math.ieee754_rem_pio2.(x))
    for _x in x
        out1 = RemPiO2.rem_pio2(_x)
        out2 = Base.Math.ieee754_rem_pio2(_x)
        @test out1[1] == out2[1]
        @test out1[2].hi == out2[2][1]
        @test out1[2].lo == out2[2][2]
    end

    @printf "%.5f    %.5f    %.5f\n" time_my_rem[i] time_rem[i] time_my_rem[i]/time_rem[i]
end
println("--"^13)
#
x = 0.32423423
xʰ = RemPiO2.highword(x)
# positive part of highword
xʰ⁺ = xʰ&0x7fffffff
RemPiO2.rem_pio2(x, xʰ, xʰ⁺)

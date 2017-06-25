using RemPiO2
using Base.Test
import Base: TwicePrecision, convert

# Welcome to the Reduction Lab

# First
x = 0.32423423
out1 = RemPiO2.rem_pio2_kernel(x)
out2 = Base.Math.ieee754_rem_pio2(x)

@test out1[1] == out2[1]
@test out1[2].hi == out2[2][1]
@test out1[2].lo == out2[2][2]

using BenchmarkTools
intervals = []
# Testing some important intervals and values around bad cases that need higher precision
for i = 1:9
    push!(intervals, linspace(-pi*i/4, -pi*(i-1)/4, 100000))
    push!(intervals, linspace(pi*(i-1)/4, pi*i/4, 100000))
end

#push!(intervals, linspace(-2.0^20, -pi*7/4, 100000))
#push!(intervals, linspace(pi*7/4, 2.0^20, 100000))
#push!(intervals, linspace(2.0^20, 2.0^60, 100000))
#intervals = [intervals[end]]
time_my_rem = []
time_rem = []
function printheader()
    println("--"^29)
    println("Testing speed and accuracy of rempio2 in [-pi*9/4, pi*9/4]")
    println("Numbers below are elapsed time returned from @belapsed")
    println()
    println("Every number is checked between the two implementations")
    println("using a @test such that any difference results in termi-")
    println("nation of the program.")
    println()
    println("(1) -> New rempio2")
    println("(2) -> Old rempio2")
    println("(3) -> (1)/(2)")
    println()
    println("interval    (1)        (2)        (3)      ")
    println("-"^41)
end

for (i, inter) in enumerate(intervals)
    x = collect(inter)
    push!(time_my_rem, @belapsed RemPiO2.rem_pio2_kernel.($x))
    push!(time_rem, @belapsed Base.Math.ieee754_rem_pio2.($x))
    for _x in x
        out1 = RemPiO2.rem_pio2_kernel(_x)
        out2 = Base.Math.ieee754_rem_pio2(_x)
        @test out1[1] == out2[1]
        @test out1[2].hi == out2[2][1]
        @test out1[2].lo == out2[2][2]
    end
    i == 1 && println("Testing values in [-pi*9/4, pi*9/4]")
    i == 1 && printheader()
    @printf "%3.0d         %.5f    %.5f    %.5f\n" i time_my_rem[i] time_rem[i] time_my_rem[i]/time_rem[i]
end



intervals = []
push!(intervals, linspace(-2.0^5, -pi*9/4, 100000))
push!(intervals, linspace(pi*9/4, 2.0^5, 100000))
push!(intervals, linspace(2.0^5, 2.0^10, 100000))
push!(intervals, linspace(-2.0^10, -2.0^5, 100000))
push!(intervals, linspace(2.0^10, 2.0^15, 100000))
push!(intervals, linspace(-2.0^15, -2.0^10, 100000))
push!(intervals, linspace(2.0^20, 2.0^15, 100000))
push!(intervals, linspace(-2.0^15, -2.0^20, 100000))

println("Testing values in [-22.0^20, -pi*9/4] and [pi*9/4, 2.0^20]")
for (i, inter) in enumerate(intervals)
    x = collect(inter)
    push!(time_my_rem, @belapsed RemPiO2.rem_pio2_kernel.($x))
    push!(time_rem, @belapsed Base.Math.ieee754_rem_pio2.($x))
    for _x in x
        out1 = RemPiO2.rem_pio2_kernel(_x)
        out2 = Base.Math.ieee754_rem_pio2(_x)
        @test out1[1] == out2[1]
        @test out1[2].hi == out2[2][1]
        @test out1[2].lo == out2[2][2]
    end
    @printf "%3.0d         %.5f    %.5f    %.5f\n" i time_my_rem[i] time_rem[i] time_my_rem[i]/time_rem[i]
end


intervals = []
for i = 1:10
    push!(intervals, linspace(2.0^(20+(i-1)*5), 2.0^(20+i*5), 100000))
    push!(intervals, linspace(-2.0^(20+(i-1)*5), -2.0^(20+i*5), 100000))
end

println("Testing values in [-2.0^70, -pi*9/4] and [pi*9/4, 2.0^70]")
for (i, inter) in enumerate(intervals)
    x = collect(inter)
    push!(time_my_rem, @belapsed RemPiO2.rem_pio2_kernel.($x))
    push!(time_rem, @belapsed Base.Math.ieee754_rem_pio2.($x))
    for _x in x
        out1 = RemPiO2.rem_pio2_kernel(_x)
        out2 = Base.Math.ieee754_rem_pio2(_x)
#        @test out1[1] == out2[1]
#        @test out1[2].hi == out2[2][1]
#        @test out1[2].lo == out2[2][2]
        @test out1[2].hi+out1[2].lo == sum(out2[2])
    end
    @printf "%3.0d         %.5f    %.5f    %.5f\n" i time_my_rem[i] time_rem[i] time_my_rem[i]/time_rem[i]
end







x = intervals[1][1]
xʰ = RemPiO2.highword(x)
# positive part of highword
xʰ⁺ = xʰ&0x7fffffff
RemPiO2.rem_pio2_kernel(x, xʰ, xʰ⁺)




# Additional notes:
# - worst case (Muller et al, 2005, p 183) is x = 6381956970095103.0 * 2.0^797
#   for which abs(rem(x,π/2, RoundNearest)) ≥ 2^-60.9
# - error of w is ~ 2^-128
# - error of f is ~ 2^-126
# - error of z is ~ |z|*2^-(53+26) + 2^-126
# - error of y is ~ 3*|y|*2^-(53+26) + 1.6*2^-126
# (TODO: check these)
# Total error ~ 2^-12 ulps

function rem_pio2(x::BigFloat)
    setprecision(BigFloat, 4096)
    rem(big(x), big(pi)/big(2), RoundNearest)
end
function big_to_hilo(x::BigFloat)
    hi = convert(Float64, x)
    lo = convert(Float64, x - hi)
    Base.TwicePrecision(hi, lo)
end

function convert(::Type{TwicePrecision{Float64}}, x::BigFloat)
    hi = convert(Float64, x)
    lo = convert(Float64, x - hi)
    Base.TwicePrecision(hi, lo)
end

x = 3.334234234
bigy = rem_pio2(big(x))
bigytp = convert(TwicePrecision{Float64}, bigy)
newy = RemPiO2.rem_pio2(x)
oldy = Base.Math.ieee754_rem_pio2(x)
oldytp = (oldy[1], TwicePrecision(oldy[2]...))

using Openlibm
Openlibm.sin(big_to_hilo(big(x)))

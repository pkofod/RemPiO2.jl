module RemPiO2

import Base: TwicePrecision, significand_mask
import Base.Math: significand_bits, exponent_mask, exponent_bias

include("codywaite.jl")
include("paynehanek.jl")

@inline highword(x::UInt64) = unsafe_trunc(UInt32,x >> 32)
@inline highword(x::Float64) = unsafe_trunc(UInt32,highword(reinterpret(UInt64, x)))

const invpio2 =  6.36619772367581382433e-01 # 0x3FE45F30, 0x6DC9C883 */
const pio2_1  =  1.57079632673412561417e+00 # 0x3FF921FB, 0x54400000 */
const pio2_1t =  6.07710050650619224932e-11 # 0x3DD0B461, 0x1A626331 */
const pio2_2  =  6.07710050630396597660e-11 # 0x3DD0B461, 0x1A600000 */
const pio2_2t =  2.02226624879595063154e-21 # 0x3BA3198A, 0x2E037073 */
const pio2_3  =  2.02226624871116645580e-21 # 0x3BA3198A, 0x2E000000 */
const pio2_3t =  8.47842766036889956997e-32 # 0x397B839A, 0x252049C1 */
#=
zero =  0.00000000000000000000e+00, /* 0x00000000, 0x00000000 */
invpio2 =  6.36619772367581382433e-01, /* 0x3FE45F30, 0x6DC9C883 */
pio2_3  =  2.02226624871116645580e-21, /* 0x3BA3198A, 0x2E000000 */
pio2_3t = 8.47842766036889956997e-32; /* 0x397B839A, 0x252049C1 */
=#


const two24 = 1.67772160000000000000e+07    # 0x41700000, 0x00000000

"""
    rem_pio2(x::Float64)

returns the remainer of x modulo π/2 as a TwicePrecision number, along with a k such
that.
"""
function rem_pio2(x::Float64)
    # get highword of x
    xʰ = highword(x)
    # positive part of highword
    xʰ⁺ = xʰ&0x7fffffff
    #	if xʰ⁺ <= 0x3fe921fb #|x| ~<= pi/4, no need for reduction, just return input
    #        return Int(0), TwicePrecision(x, nothing)
    #    end
    _rem_pio2_kernel(x, xʰ, xʰ⁺)
end

#"""
#    _rem_pio2_kernel(x, xʰ, xʰ⁺)

#returns the remainder of x modulo π/2 as a TwicePrecision number, along with a k
#such that.
#"""

function _rem_pio2_kernel(x, xʰ, xʰ⁺)
    if xʰ⁺ <= 0x400f6a7a #|x| ~<= 5pi/4, use Cody Waite with two constants
        # |x| ~= pi/2 or 2pi/2, use precise Cody Waite schemee
        if (xʰ⁺ & 0xfffff) == 0x921fb

            return cody_waite_3c(x, xʰ⁺)
        end
        # |x| ~<= 3pi/4
        if xʰ⁺ <= 0x4002d97c
            if x > 0
                return 1, cody_waite_2c(x, 1.0)
            else
                return -1, cody_waite_2c(x, -1.0)
            end
        else
            if x > 0
                return 2, cody_waite_2c(x, 2.0)
            else
                return -2, cody_waite_2c(x, -2.0)
            end
        end
    end

    if xʰ⁺ <= 0x401c463b # |x| ~<= 9pi/4, use Cody Waite with two constants
        if (xʰ⁺ <= 0x4015fdbc) # |x| ~<= 7pi/4
            if (xʰ⁺ == 0x4012d97c) # |x| ~= 3pi/2, use precise Cody Waite scheme
                return cody_waite_3c(x, xʰ⁺)
            end
            if x > 0
                return 3, cody_waite_2c(x, 3.0)
            else
                return -3, cody_waite_2c(x, -3.0)
            end
        else
            if xʰ⁺ == 0x401921fb # |x| ~= 4pi/2, use precise Cody Waite scheme
                return cody_waite_3c(x, xʰ⁺)
            end
            if x > 0
                return 4, cody_waite_2c(x, 4.0)
            else
                return -4, cody_waite_2c(x, -4.0)
            end
        end
    end
    if xʰ⁺<0x413921fb # |x| ~< 2^20*pi/2, use precise Cody Waite scheme
        return cody_waite_3c(x, xʰ⁺)
    end

    # if |x| >= 2^20*pi/2 switch to Payne Hanek
    return paynehanek(x)
end

end # module

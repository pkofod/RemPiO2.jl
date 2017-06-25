module RemPiO2

import Base: TwicePrecision, significand_mask
import Base.Math: significand_bits, exponent_mask, exponent_bias

include("codywaite.jl")
include("paynehanek.jl")

"""
    highword(x)

returns the high word of x as a UInt32.
"""
@inline highword(x::UInt64) = unsafe_trunc(UInt32,x >> 32)
@inline highword(x::Float64) = unsafe_trunc(UInt32,highword(reinterpret(UInt64, x)))

"""
    rem_pio2(x::Float64)

returns the remainder of x modulo π/2 as a TwicePrecision number, along with a k
the specifies
"""

function rem_pio2(x::Float64)
    xh = highword(x)
    xhp = xh & 0x7fffffff # positive part of highword

    if xhp <= 0x3fe921fb #|x| ~<= pi/4, no need for reduction, just return input
        return Int(0), TwicePrecision(x, 0.0)
    end
    rem_pio2_kernel(x, xh, xhp)
end

"""
    rem_pio2_kernel(x, xh, xhp)

returns the remainder of x modulo π/2 as a TwicePrecision number, along with a k
such that.
"""
function rem_pio2_kernel(x)
    xh = highword(x)
    xhp = xh & 0x7fffffff # positive part of highword
    rem_pio2_kernel(x, xh, xhp)
end
function rem_pio2_kernel(x, xh, xhp)
    if xhp <= 0x400f6a7a #|x| ~<= 5pi/4, use Cody Waite with two constants
        if (xhp & 0xfffff) == 0x921fb # |x| ~= pi/2 or 2pi/2, use precise Cody Waite scheme
            return cody_waite_ext_pio2(x, xhp)
        end
        if xhp <= 0x4002d97c # |x| ~<= 3pi/4
            if x > 0
                return 1, cody_waite_2c_pio2(x, 1.0)
            else
                return -1, cody_waite_2c_pio2(x, -1.0)
            end
        else
            if x > 0
                return 2, cody_waite_2c_pio2(x, 2.0)
            else
                return -2, cody_waite_2c_pio2(x, -2.0)
            end
        end
    end

    if xhp <= 0x401c463b # |x| ~<= 9pi/4, use Cody Waite with two constants
        if (xhp <= 0x4015fdbc) # |x| ~<= 7pi/4
            if (xhp == 0x4012d97c) # |x| ~= 3pi/2, use precise Cody Waite scheme
                return cody_waite_ext_pio2(x, xhp)
            end
            if x > 0
                return 3, cody_waite_2c_pio2(x, 3.0)
            else
                return -3, cody_waite_2c_pio2(x, -3.0)
            end
        else
            if xhp == 0x401921fb # |x| ~= 4pi/2, use precise Cody Waite scheme
                return cody_waite_ext_pio2(x, xhp)
            end
            if x > 0
                return 4, cody_waite_2c_pio2(x, 4.0)
            else
                return -4, cody_waite_2c_pio2(x, -4.0)
            end
        end
    end
    if xhp<0x413921fb # |x| ~< 2^20*pi/2, use precise Cody Waite scheme
        return cody_waite_ext_pio2(x, xhp)
    end

    # if |x| >= 2^20*pi/2 switch to Payne Hanek
    return paynehanek(x)
end

end # module

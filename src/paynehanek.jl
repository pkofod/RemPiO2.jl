const INV2PI = UInt64[
    0x28be_60db_9391_054a,
    0x7f09_d5f4_7d4d_3770,
    0x36d8_a566_4f10_e410,
    0x7f94_58ea_f7ae_f158,
    0x6dc9_1b8e_9093_74b8,
    0x0192_4bba_8274_6487,
    0x3f87_7ac7_2c4a_69cf,
    0xba20_8d7d_4bae_d121,
    0x3a67_1c09_ad17_df90,
    0x4e64_758e_60d4_ce7d,
    0x2721_17e2_ef7e_4a0e,
    0xc7fe_25ff_f781_6603,
    0xfbcb_c462_d682_9b47,
    0xdb4d_9fb3_c9f2_c26d,
    0xd3d1_8fd9_a797_fa8b,
    0x5d49_eeb1_faf9_7c5e,
    0xcf41_ce7d_e294_a4ba,
    0x9afe_d7ec_47e3_5742,
    0x1580_cc11_bf1e_daea]

function paynehanek(x::Float64)
    # 1. Convert to form
    #
    #    x = X * 2^k,
    #
    # where 2^(n-1) <= X < 2^n  is an n-bit integer (n = 53, k = exponent(x)-52 )

    u = reinterpret(UInt64, x)
    X = (u & significand_mask(Float64)) | (one(UInt64) << significand_bits(Float64))
    k = Int((u & exponent_mask(Float64)) >> significand_bits(Float64)) - exponent_bias(Float64) - significand_bits(Float64)

    # 2. Let α = 1/2π, then:
    #
    #    α*x mod 1 ≡ [(α*2^k mod 1)*X] mod 1
    #
    # so we can ignore the first k bits of α. Extract the next 3 64-bit parts of α.
    #
    # i.e. equivalent to
    #     setprecision(BigFloat,4096)
    #     α  = 1/(2*big(pi))
    #     A  = mod(ldexp(α,k), 1)
    #     z1 = ldexp(A,64)
    #     a1 = trunc(UInt64, z1)
    #     z2 = ldexp(z1-a1, 64)
    #     a2 = trunc(UInt64, z2)
    #     z3 = ldexp(z2-a2, 64)
    #     a3 = trunc(UInt64, z3)

    # idx, shift = divrem(k, 64), but divrem is slower
    idx = k >> 6
    shift = k - (idx << 6)
    if shift == 0
        a1 = INV2PI[idx+1]
        a2 = INV2PI[idx+2]
        a3 = INV2PI[idx+3]
    else
        # use shifts to extract the relevant 64 bit window
        a1 = (idx < 0 ? zero(UInt64) : INV2PI[idx+1] << shift) | (INV2PI[idx+2] >> (64 - shift))
        a2 = (INV2PI[idx+2] << shift) | (INV2PI[idx+3] >> (64 - shift))
        a3 = (INV2PI[idx+3] << shift) | (INV2PI[idx+4] >> (64 - shift))
    end

    # 3. Perform the multiplication:
    #
    #      X.  0  0  0
    #   ×  0. a1 a2 a3
    #   ==============
    #      _.  w  w  _
    #
    # (i.e. ignoring integer and lowest bit parts of result)

    w1 = UInt128(X*a1) << 64 # overflow becomes integer
    w2 = widemul(X,a2)
    w3 = widemul(X,a3) >> 64
    w = w1 + w2 + w3         # quotient fraction after division by 2π

    # adjust for sign of x
    w = flipsign(w,x)

    # 4. convert to quadrant, quotient fraction after division by π/2:
    q = (((w>>125)%Int +1)>>1) # nearest quadrant
    f = (w<<2) % Int128 # fraction part of quotient after division by π/2, taking values on [-0.5,0.5)

    # 5. convert quotient fraction to split precision Float64
    z_hi,z_lo = fromfraction(f)

    # 6. multiply by π/2
    pio2_hi = 1.5707963407039642
    pio2_lo = -1.3909067614167116e-8

    y_hi = (z_hi+z_lo)*(pio2_hi+pio2_lo)
    y_lo = (((z_hi*pio2_hi - y_hi) + z_hi*pio2_lo) + z_lo*pio2_hi) + z_lo*pio2_lo
    return q, TwicePrecision(y_hi, y_lo)
end

"""
    fromfraction(f::Int128)
Computes a tuple of values `(y1,y2)` such that
    y1 + y2 == f / 2^128
and the significand of `y1` has 27 trailing zeros.
"""
function fromfraction(f::Int128)
    if f == 0
        return (0.0,0.0)
    end

    # 1. get leading term truncated to 26 bits
    s = ((f < 0) % UInt64) << 63     # sign bit
    x = abs(f) % UInt128             # magnitude
    n1 = 128-leading_zeros(x)         # ndigits0z(x,2)
    m1 = ((x >> (n1-26)) % UInt64) << 27
    d1 = ((n1-128+1021) % UInt64) << 52
    z1 = reinterpret(Float64, s | (d1 + m1))

    # 2. compute remaining term
    x2 = (x - (UInt128(m1) << (n1-53)))
    if x2 == 0
        return (z1, 0.0)
    end
    n2 = 128-leading_zeros(x2)
    m2 = (x2 >> (n2-53)) % UInt64
    d2 = ((n2-128+1021) % UInt64) << 52
    z2 = reinterpret(Float64,  s | (d2 + m2))
    return (z1,z2)
end

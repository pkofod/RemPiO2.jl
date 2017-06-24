function cody_waite_2c(x, signed_k)
    z = x - signed_k*pio2_1;
    y1 = z - signed_k*pio2_1t;
    y2 = (z - y1) - signed_k*pio2_1t;
    TwicePrecision(y1, y2)
end

function cody_waite_3c(x, xʰ⁺)
    fn = x*invpio2+0x1.8p52
    fn = fn-0x1.8p52
    n  = Int32(fn)
    r  = x-fn*pio2_1
    w  = fn*pio2_1t # 1st round good to 85 bit
    j  = xʰ⁺>>20
    y1 = r-w
    high = highword(y1)
    i = j-((high>>20)&0x7ff)
    if i>16  # 2nd iteration needed, good to 118
        t  = r
        w  = fn*pio2_2
        r  = t-w
        w  = fn*pio2_2t-((t-r)-w)
        y1 = r-w
        high = highword(y1)
        i = j-((high>>20)&0x7ff)
        if i>49 # 3rd iteration need, 151 bits acc
            t  = r # will cover all possible cases
            w  = fn*pio2_3
            r  = t-w
            w  = fn*pio2_3t-((t-r)-w)
            y1 = r-w
        end
    end
    y2 = (r-y1)-w
    return Int(n), TwicePrecision(y1, y2)
end

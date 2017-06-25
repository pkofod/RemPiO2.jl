const invpio2 =  6.36619772367581382433e-01 # 0x3FE45F30, 0x6DC9C883 */
const pio2_1  =  1.57079632673412561417e+00 # 0x3FF921FB, 0x54400000 */
const pio2_1t =  6.07710050650619224932e-11 # 0x3DD0B461, 0x1A626331 */
const pio2_2  =  6.07710050630396597660e-11 # 0x3DD0B461, 0x1A600000 */
const pio2_2t =  2.02226624879595063154e-21 # 0x3BA3198A, 0x2E037073 */
const pio2_3  =  2.02226624871116645580e-21 # 0x3BA3198A, 0x2E000000 */
const pio2_3t =  8.47842766036889956997e-32 # 0x397B839A, 0x252049C1 */

function cody_waite_2c_pio2(x, signed_k)
    z = x - signed_k*pio2_1;
    y1 = z - signed_k*pio2_1t;
    y2 = (z - y1) - signed_k*pio2_1t;
    TwicePrecision(y1, y2)
end

function cody_waite_ext_pio2(x, xʰ⁺)
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

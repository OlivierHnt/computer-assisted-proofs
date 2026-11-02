# Reading the frozen inputs: raw little-endian IEEE binary64, no header. The
# approximate eigenvalue λ̄ is one real number (cases II, III) or the pair
# (Re λ̄, Im λ̄) (case I); it is turned into a thin interval.

function read_f64(path, nvalues)
    filesize(path) == 8nvalues || error("$(basename(path)) has the wrong size: expected $(8nvalues) bytes")
    return ltoh.(reinterpret(Float64, read(path)))
end

read_real(path) = (v = read_f64(path, 1); interval(v[1]))
read_complex(path) = (v = read_f64(path, 2); complex(interval(v[1]), interval(v[2])))

# printing: endpoints in full (shortest round-trip decimal of each binary64)
enc(x::Interval) = string("[", inf(x), ", ", sup(x), "]")
enc(z::Complex{<:Interval}) = string(enc(real(z)), " + im*", enc(imag(z)))

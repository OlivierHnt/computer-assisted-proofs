@enum ReturnCode begin
    Success
    UnboundedY      # a component of Y is not a bounded interval
    UnboundedZ      # an entry of Z (or W) is not a bounded interval
    NoRadius        # no radii with nonpositive radii polynomials and a contracting Z
    NotGuaranteed   # a bound is not guaranteed (decoration)
end

struct Certificate{P,Alg,TX,TY,TZ,TW,TR,Tr,Te}
    problem   :: P
    algorithm :: Alg
    x_bar     :: TX
    Y         :: TY   # one bound per component
    Z         :: TZ   # matrix of block norms
    W         :: TW   # second-order tensor, zero for order 1
    R         :: TR   # thresholds of the radii
    radii     :: Tr
    eta       :: Te   # test vector of the contraction
    retcode   :: ReturnCode
end

function Certificate(pb, alg, x_bar, Y, Z, W, R)
    radii, eta, success = set_of_radii(Y, Z, W, R)
    guaranteed = all(isguaranteed, Y) && all(isguaranteed, Z) && all(isguaranteed, W) && all(isguaranteed, radii)
    retcode = success & guaranteed ? Success :
              !all(_isfinite, Y) ? UnboundedY :
              !(all(_isfinite, Z) && all(_isfinite, W)) ? UnboundedZ :
              !success ? NoRadius : NotGuaranteed
    return Certificate(pb, alg, x_bar, Y, Z, W, R, radii, eta, retcode)
end
_isfinite(x) = isbounded(x) && !isempty_interval(x)

LinearAlgebra.issuccess(cert::Certificate) = cert.retcode == Success

function radius(cert::Certificate)
    issuccess(cert) || throw(ArgumentError("no radius: the certificate failed with $(cert.retcode)"))
    return inf.(cert.radii)
end

# the enclosure of the zero, componentwise, to be used as input of the next problem
function solution(cert::Certificate)
    issuccess(cert) || throw(ArgumentError("no solution: the certificate failed with $(cert.retcode)"))
    return ball(cert.problem, cert.x_bar, radius(cert))
end

function Base.show(io::IO, cert::Certificate)
    print(io, "Certificate(", cert.retcode, "; radii = ", issuccess(cert) ? radius(cert) : "none", ", Y = ", sup.(cert.Y), ")")
end

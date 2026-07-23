import SatelliteToolboxTransformations

import LinearAlgebra.cross
function cross(x::AbstractVector)
    return [0 -x[3] x[2];
        x[3] 0 -x[1];
        -x[2] x[1] 0]
end
function uncross(X::AbstractMatrix)
    return [-X[2,3]; X[1,3]; -X[1,2]]
end
function residualSO3(X::AbstractMatrix)
    if any(isnan.(X[:]))
        return Inf
    else
        return (det(X) - 1) + (sum(X'*X .- diagm([1,1,1])))
    end
end
function residualso3(X::AbstractMatrix)
    return sum(X' + X) 
end
function projSO3(X::AbstractMatrix)
    (U, S, V) = svd(X)

    C_BA = U*diagm([1, 1, det(U)*det(V)])*V'
    return C_BA
end

@doc raw"""
    r_min_arc(x_A::AbstractVector, x_B::AbstractVector)

Return a rotation matrix `C` such that `C  * x_A == x_B`.
"""
function r_min_arc(x_A::AbstractVector{<:Real}, x_B::AbstractVector{<:Real})::AbstractMatrix{<:Real}
    tol = 1e-12
    x_A = normalize(x_A)
    x_B = normalize(x_B)
    ax = zeros(3)
    cx = LinearAlgebra.cross(x_A, x_B)
    if norm(cx) > tol
        ax = normalize(cx)
    end
    cang = clamp(dot(x_A, x_B), -1, 1)
    return I(3)*cang + (1 - cang)*ax*ax' + cross(ax)*sqrt(1 - cang^2)
end

function worm(k, w)
    s = 'o'^(w - 1)
    out = s[1:k]*'O'*s[k+1:end]
    return out
end


global const eop = SatelliteToolboxTransformations.fetch_iers_eop()
import SatelliteToolboxTransformations
import Makie
import LinearAlgebra.cross

function cross(x::AbstractVector)
    return [
        0 -x[3] x[2];
        x[3] 0 -x[1];
        -x[2] x[1] 0
    ]
end
function uncross(X::AbstractMatrix)
    return [-X[2, 3]; X[1, 3]; -X[1, 2]]
end
function residualSO3(X::AbstractMatrix)
    if any(isnan.(X[:]))
        return Inf
    else
        return (det(X) - 1) + (sum(X'*X .- diagm([1, 1, 1])))
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
function r_min_arc(
    x_A::AbstractVector{<:Real},
    x_B::AbstractVector{<:Real},
)::AbstractMatrix{<:Real}
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

function dcm_to_quat_tr(X::AbstractMatrix)
    qr = 0.5 * sqrt(1.0 + LinearAlgebra.tr(X))
    qi = 1 / 4 / qr * (X[2, 3] - X[3, 2])
    qj = 1 / 4 / qr * (X[3, 1] - X[1, 3])
    qk = 1 / 4 / qr * (X[1, 2] - X[2, 1])
    return Makie.Quaternion(qi, qj, qk, qr)
end
function dcm_to_quat_11(X::AbstractMatrix)
    qi = 0.5 * sqrt(1.0 + X[1, 1] - X[2, 2] - X[3, 3])
    qr = 1 / 4 / qi * (X[2, 3] - X[3, 2])
    qj = 1 / 4 / qi * (X[1, 2] + X[2, 1])
    qk = 1 / 4 / qi * (X[3, 1] + X[1, 3])
    return Makie.Quaternion(qi, qj, qk, qr)
end
function dcm_to_quat_22(X::AbstractMatrix)
    qj = 0.5 * sqrt(1.0 - X[1, 1] + X[2, 2] - X[3, 3])
    qr = 1 / 4 / qj * (X[3, 1] - X[1, 3])
    qi = 1 / 4 / qj * (X[1, 2] + X[2, 1])
    qk = 1 / 4 / qj * (X[2, 3] + X[3, 2])
    return Makie.Quaternion(qi, qj, qk, qr)
end
function dcm_to_quat_33(X::AbstractMatrix)
    qk = 0.5 * sqrt(1.0 - X[1, 1] - X[2, 2] + X[3, 3])
    qr = 1 / 4 / qk * (X[1, 2] - X[2, 1])
    qi = 1 / 4 / qk * (X[3, 1] + X[1, 3])
    qj = 1 / 4 / qk * (X[2, 3] + X[3, 2])
    return Makie.Quaternion(qi, qj, qk, qr)
end
function dcm_to_quat(X::AbstractMatrix)
    # implements Shepperd 1978 method for choosing a numerically favorable conversion.
    (v, p) = findmax([LinearAlgebra.tr(X), X[1, 1], X[2, 2], X[3, 3]])
    if p == 1
        return dcm_to_quat_tr(X)
    elseif p == 2
        return dcm_to_quat_11(X)
    elseif p == 3
        return dcm_to_quat_22(X)
    elseif p == 4
        return dcm_to_quat_33(X)
    end
end
Makie.Quaternion(X::AbstractMatrix) = begin
    qr = 0.5 * sqrt(1 + LinearAlgebra.tr(X))
    qi = 1 / 4 / qr * (X[3, 2] - X[2, 3])
    qj = 1 / 4 / qr * (X[1, 3] - X[3, 1])
    qk = 1 / 4 / qr * (X[2, 1] - X[1, 2])
    return Makie.Quaternion(qi, qj, qk, qr)
end

@doc raw"""
    worm(k::Int, w::Int)

Makea little string of width `w`, and accent the `k`th character.

Use to print a spinner in a loop by passing some some counter `c % w` as `k`.
"""
function worm(k::Int, w::Int)
    s = 'o'^(w - 1)
    out = s[1:k]*'O'*s[(k+1):end]
    return out
end

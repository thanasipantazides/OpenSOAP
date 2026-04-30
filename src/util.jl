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
    return (det(X) - 1) + (sum(X'*X .- diagm([1,1,1])))
end
function residualso3(X::AbstractMatrix)
    return sum(X' + X) 
end
function projSO3(X::AbstractMatrix)
    (U, S, V) = svd(X)

    C_BA = U*diagm([1, 1, det(U)*det(V)])*V'
    return C_BA
end
function r_min_arc(x_A::AbstractVector{<:Real}, x_B::AbstractVector{<:Real})::AbstractMatrix{<:Real}
    ax = normalize(LinearAlgebra.cross(x_A, x_B))
    cang = dot(x_A, x_B) / norm(x_A) / norm(x_B)
    return I(3)*cang + (1 - cang)*ax*ax' + cross(ax)*sqrt(1 - cang^2)
end

function worm(k, w)
    s = 'o'^(w - 1)
    out = s[1:k]*'O'*s[k+1:end]
    return out
end


global const eop = SatelliteToolboxTransformations.fetch_iers_eop()
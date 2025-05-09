using StaticArrays
import LinearAlgebra, LinearAlgebra.cross
@doc raw"""
    cross(x)

Return the skew-symmetric matrix implementing the `cross` (``\times``) operator for the provided 3-vector `x`. For a 3-vectors `x` and `y`,
`cross(x)*y == LinearAlgebra.cross(x, y)`

This is useful, e.g., when working in the Lie group ``\mathfrak{so}(3)``, the tangent space to the rotation group ``\mathrm{SO}(3)``.
"""
function LinearAlgebra.cross(x::Vector{<:Number})::Matrix{<:Number}
    return [
        0 -x[3] x[2];
        x[3] 0 -x[1];
        -x[2] x[1] 0
    ]
end

function LinearAlgebra.cross(x::SVector{3,Float64})::SMatrix{3,3,Float64}
    return [
        0 -x[3] x[2];
        x[3] 0 -x[1];
        -x[2] x[1] 0
    ]
end

function axis(X::SMatrix{3,3})::SVector{3}
    return uncross(X - X')
end
function axis(X::Matrix{<:Real})::Vector{<:Real}
    return uncross(X - X')
end
function axisangle(X::Matrix{<:Real})
    λ, V = LinearAlgebra.eigen(Matrix(X))
    i = sortperm(λ, by=imag)
    ax = real.(V[:, i[2]])
    # ang = atan(imag(λ[i[3]]), real(λ[i[3]]))

    ang = acos((tr(X) - 1) / 2)
    e1 = (X[3, 2] - X[2, 3]) / 2 / sin(ang)
    e2 = (X[1, 3] - X[3, 1]) / 2 / sin(ang)
    e3 = (X[2, 1] - X[1, 2]) / 2 / sin(ang)
    ax = [e1, e2, e3]
    return (ax, ang)


end
function axisangle(X::SMatrix{3,3})
    skew = X - X'
    ax = uncross(skew)
    λ, V = LinearAlgebra.eigen(Matrix(X))
    i = sortperm(λ, by=imag)
    ax = real.(V[:, i[2]])
    ang = asin(min(1.0, -imag(λ[i[1]])))
    return (ax, ang)
end
function axisangle(ax::Vector{<:Real}, ang::Real)
    return I * cos(ang) + (1 - cos(ang)) * ax * ax' + cross(ax) * sin(ang)
end

@doc raw"""
    uncross(X)

Return the vector `u` that has been turned into a matrix cross operator `X`, such that ``Xv = \texttt{uncross}(X)\times v``.

This is useful, e.g., when working in the Lie group ``\mathfrak{so}(3)``, the tangent space to the rotation group ``\mathrm{SO}(3)``.
"""
function uncross(X::Matrix{<:Number})::Vector{<:Number}
    u = [-X[2, 3]; X[1, 3]; -X[1, 2]]
    if all(u .== 0)
        return zeros(3)
    else
        return u / norm(u)
    end
end
function uncross(X::SMatrix{3,3})::SVector{3}
    u = [-X[2, 3]; X[1, 3]; -X[1, 2]]
    if all(u .== 0)
        return zeros(3)
    else
        return u / norm(u)
    end
    return u / norm(u)
end


@doc raw"""
    r_min_arc(x_A, x_B)

Return the DCM implementing the minimum arc-length rotation that maps vector `x_A` into vector `x_B`. This uses the axis-angle parameterization to compute the DCM.
"""
function r_min_arc(x_A::Vector{<:Real}, x_B::Vector{<:Real})::Matrix{<:Real}
    ax = LinearAlgebra.cross(x_A, x_B) / norm(x_A) / norm(x_B)
    ax = ax / norm(ax)
    cang = x_A' * x_B / norm(x_A) / norm(x_B)

    return I * cang + (1 - cang) * ax * ax' + cross(ax) * sqrt(1 - cang^2)
end

function r_random()::Matrix{<:Real}
    ax = 2 * rand(3) .- 1
    ax = ax / norm(ax)
    cang = 2 * rand() - 1

    return I * cang + (1 - cang) * ax * ax' + cross(ax) * sqrt(1 - cang^2)
end

function r_euler3(ang::Real)::Matrix{<:Real}
    return [cos(ang) -sin(ang) 0;
        sin(ang) cos(ang) 0;
        0 0 1]
end
function r_euler2(ang::Real)::Matrix{<:Real}
    return [cos(ang) 0 sin(ang);
        -sin(ang) 0 cos(ang);
        0 1 0]
end
function r_euler1(ang::Real)::Matrix{<:Real}
    return [1 0 0;
        0 cos(ang) -sin(ang);
        0 sin(ang) cos(ang)]
end
function ang321(C::Matrix{<:Real})::Vector{<:Real}
    # z - y - X
    # yaw - pitch - roll
    # pitch = -asin(C[3,1])
    # cosp = cos(pitch)
    # roll = asin(C[3,2] / cosp)
    # yaw = acos(C[1,1] / cosp)
    yaw = atan(C[1, 2] / C[1, 1])
    pitch = -asin(C[1, 3])
    roll = atan(C[2, 3] / C[3, 3])
    return [yaw; pitch; roll]
end

function ang321(C::SMatrix{3,3})::SVector{3}
    # z - y - X
    # yaw - pitch - roll
    yaw = atan(C[1, 2] / C[1, 1])
    pitch = -asin(C[1, 3])
    roll = atan(C[2, 3] / C[3, 3])
    return [yaw; pitch; roll]
end

@doc raw"""
  rotinterp(R0, Rf, n)

Interpolate rotation matrices between ``R0`` and ``Rf``, in ``n`` steps.
"""
function rotinterp(R0, Rf, n)
    # Rf0 = (Rf' * R0)
    Rf0 = R0' * Rf
    (ax, ang) = axisangle(Rf0)
    # println(ax)
    # println(ang)

    R = Array{typeof(R0[1]),3}(undef, 3, 3, n)
    for k in 1:n
        kang = ang * (k - 1) / n
        R[:, :, k] = axisangle(ax, kang)
    end

    return R
end

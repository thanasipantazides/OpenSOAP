import LinearAlgebra
@doc raw"""
    cross(x)

Return the skew-symmetric matrix implementing the `cross` (``\times``) operator for the provided 3-vector `x`. For a 3-vectors `x` and `y`,
`cross(x)*y == LinearAlgebra.cross(x, y)`

This is useful, e.g., when working in the Lie group ``\mathfrak{so}(3)``, the tangent space to the rotation group ``\mathrm{SO}(3)``.
"""
function cross(x::Vector{<:Number})::Matrix{<:Number}
    return [
        0   -x[3]    x[2];
        x[3]    0   -x[1];
        -x[2]   x[1]    0;
    ]
end

@doc raw"""
    r_min_arc(x_A, x_B)

Return the DCM implementing the minimum arc-length rotation that maps vector `x_A` into vector `x_B`. This uses the axis-angle parameterization to compute the DCM.
"""
function r_min_arc(x_A::Vector{<:Real}, x_B::Vector{<:Real})::Matrix{<:Real}
    ax = LinearAlgebra.cross(x_A, x_B)/norm(x_A)/norm(x_B)
    ax = ax/norm(ax)
    cang = x_A'*x_B/norm(x_A)/norm(x_B)

    return I*cang + (1 - cang)*ax*ax' + cross(ax)*sqrt(1 - cang^2)
end

function r_random()::Matrix{<:Real}
    ax = rand(3)
    ax = ax/norm(ax)
    cang = 2*rand() - 1

    return I*cang + (1 - cang)*ax*ax' + cross(ax)*sqrt(1 - cang^2)
end
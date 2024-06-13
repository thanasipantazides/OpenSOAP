import LinearAlgebra
function cross(x::Vector{<:Number})::Matrix{<:Number}
    return [
        0   -x[3]    x[2];
        x[3]    0   -x[1];
        -x[2]   x[1]    0;
    ]
end

function r_min_arc(x_A::Vector{<:Real}, x_B::Vector{<:Real})::Matrix{<:Real}
    ax = LinearAlgebra.cross(x_A, x_B)/norm(x_A)/norm(x_B)
    cang = x_A'*x_B/norm(x_A)/norm(x_B)

    return I*cang + (1 - cang)ax*ax' + cross(ax)*sqrt(1 - cang^2)
end
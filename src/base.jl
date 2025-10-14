using StaticArrays
import LinearAlgebra, LinearAlgebra.cross
import Makie
import Base: +, *

# -------------------------------------------------------------------
# Spacecraft operational modes (used throughout simulation)
# -------------------------------------------------------------------
@enum Modes begin 
    safe
    idle
    detumble
    pointing
    charging
    science
    downlink
end

# -------------------------------------------------------------------
# Spacecraft state (dynamic + system parameters)
# -------------------------------------------------------------------
mutable struct State{T<:Real}
    position::SVector{3,T}         # spacecraft position (m)
    velocity::SVector{3,T}         # spacecraft velocity (m/s)
    angular_velocity::SVector{3,T} # spacecraft angular velocity (rad/s)
    attitude::SMatrix{3,3,T}       # direction cosine matrix (body→inertial)
    battery::T                     # current battery energy (J or Wh)
    storage::T                     # onboard data storage (bytes)
    mode::Modes                    # operational mode (enum type)
end

# -------------------------------------------------------------------
# Constructors
# -------------------------------------------------------------------

# Create a State from StaticArrays inputs
State{S}(pos::SVector{3,S}, vel::SVector{3,S}, ang_vel::SVector{3,S},
          att::SMatrix{3,3,S}, batt::S, stor::S,
          mod::Union{S,Int64}) where {S<:Real} =
    State{S}(pos, vel, ang_vel, att, batt, stor, Modes(Int(round(mod))))

# Create a State from standard Arrays (converts them to StaticArrays)
State{S}(pos::Vector{S}, vel::Vector{S}, ang_vel::Vector{S},
          att::Matrix{S}, batt::S, stor::S,
          mod::Union{S,Int64}) where {S<:Real} =
    State{S}(SVector{3}(pos[1:3]),
              SVector{3}(vel[1:3]),
              SVector{3}(ang_vel[1:3]),
              SMatrix{3,3}(att[1:3,1:3]),
              batt, stor, Modes(Int(round(mod))))

# Empty/default State
State{S}() where {S<:Real} =
    State{S}(zeros(S,3), zeros(S,3), zeros(S,3),
              SMatrix{3,3}(I), zero(S), zero(S), safe)

# -------------------------------------------------------------------
# Arithmetic Operators
# -------------------------------------------------------------------
function +(a::State, b::State)
    return State{typeof(a.position[1])}(
        a.position .+ b.position,
        a.velocity .+ b.velocity,
        a.angular_velocity .+ b.angular_velocity,
        a.attitude .+ b.attitude,
        a.battery .+ b.battery,
        a.storage .+ b.storage,
        a.mode
    )
end

function *(a::Real, b::State)
    return State{typeof(b.position[1])}(
        a .* b.position,
        a .* b.velocity,
        a .* b.angular_velocity,
        a .* b.attitude,
        a * b.battery,
        a * b.storage,
        b.mode
    )
end

# -------------------------------------------------------------------
# Cross product and quaternion utilities
# -------------------------------------------------------------------
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
    ang = acos(min(1.0, (tr(X) - 1) / 2))
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

Makie.Quaternion(X::Matrix{Float64}) = begin
    qr = 0.5 * sqrt(1 + LinearAlgebra.tr(X))
    qi = 1 / 4 / qr * (X[3, 2] - X[2, 3])
    qj = 1 / 4 / qr * (X[1, 3] - X[3, 1])
    qk = 1 / 4 / qr * (X[2, 1] - X[1, 2])
    return Makie.Quaternion(qi, qj, qk, qr)
end

Makie.Quaternion(X::SMatrix{3,3}) = begin 
    qr = 0.5 * sqrt(1 + LinearAlgebra.tr(X))
    qi = 1 / 4 / qr * (X[3, 2] - X[2, 3])
    qj = 1 / 4 / qr * (X[1, 3] - X[3, 1])
    qk = 1 / 4 / qr * (X[2, 1] - X[1, 2])
    return Makie.Quaternion(qi, qj, qk, qr)
end

# -------------------------------------------------------------------
# Uncross and rotation helpers
# -------------------------------------------------------------------
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
end

function r_min_arc(x_A::AbstractVector{<:Real}, x_B::AbstractVector{<:Real})::Matrix{<:Real}
    ax = LinearAlgebra.cross(x_A, x_B)
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

# Euler rotations
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

# Euler 3-2-1 (yaw-pitch-roll)
function ang321(C::Matrix{<:Real})::Vector{<:Real}
    yaw = atan(C[1, 2] / C[1, 1])
    pitch = -asin(C[1, 3])
    roll = atan(C[2, 3] / C[3, 3])
    return [yaw; pitch; roll]
end

function ang321(C::SMatrix{3,3})::SVector{3}
    yaw = atan(C[1, 2] / C[1, 1])
    pitch = -asin(C[1, 3])
    roll = atan(C[2, 3] / C[3, 3])
    return [yaw; pitch; roll]
end

# -------------------------------------------------------------------
# Interpolation functions
# -------------------------------------------------------------------
function rotinterp(R0::AbstractMatrix{<:Real}, Rf::AbstractMatrix{<:Real}, n::Int)
    Rf0 = R0' * Rf
    (ax, ang) = axisangle(Rf0)
    R = Array{typeof(R0[1]),3}(undef, 3, 3, n)
    for k in 1:n
        kang = ang * (k - 1) / n
        R[:, :, k] = axisangle(ax, kang)
    end
    return R
end

function lininterp(v0::AbstractVector{<:Real}, vf::AbstractVector{<:Real}, n::Int)
    m = length(v0)
    vdir = vf - v0
    V = Matrix{typeof(v0[1])}(undef, m, n)
    for k in 1:n
        V[:,k] = v0 + (k - 1) / n * vdir
    end
    return V
end

function interp(p0::Real, pf::Real, n::Int)
    return [p0 + (pf - p0)*(k - 1)/n for k in 1:n]
end

function interp(s0::State{<:Real}, sf::State{<:Real}, n::Int)
    T = typeof(s0.position[1])
    s = Vector{State{T}}()

    positions = lininterp(s0.position, sf.position, n)
    velocities = lininterp(s0.velocity, sf.velocity, n)
    angular_velocities = lininterp(s0.angular_velocity, sf.angular_velocity, n)
    attitudes = rotinterp(s0.attitude, sf.attitude, n)
    batteries = interp(s0.battery, sf.battery, n)
    storages = interp(s0.storage, sf.storage, n)
    modes = fill(s0.mode, n)

    for k in 1:n
        push!(s, State{T}(positions[:,k], velocities[:,k], angular_velocities[:,k],
                          attitudes[:,:,k], batteries[k], storages[k], modes[k]))
    end
    return s
end

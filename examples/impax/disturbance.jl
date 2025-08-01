# === Drag Force and Torque ===
# This function estimates the drag force and the resulting torque on the spacecraft
# due to atmospheric drag.
function compute_drag(r, v)
    ad = 1e-12                # Atmospheric density at orbit altitude (kg/m^3)
    Cd = 2.2                  # Drag coefficient... how "aerodynamic" the shape is
    A = 1.0                   # Cross-sectional area (m^2) of the spacecraft facing the velocity direction
    v_rel = v                 # Velocity of the spacecraft relative to the atmosphere

    # Compute the drag force using the drag equation:
    # F_drag = -0.5 * density * velocity^2 * drag coefficient * area
    F_drag = -0.5 * ad * norm(v_rel)^2 * Cd * A * normalize(v_rel)

    # Compute the torque caused by drag force acting away from the center of mass
    T_drag = cross([0.1, 0.0, 0.0], F_drag)

    # Return both the drag force and the resulting torque
    return F_drag, T_drag
end

# === Solar Radiation Pressure ===
# This function estimates the force and torque from sunlight pushing on the spacecraft.
function compute_solar_pressure()
    P = 4.5e-6                # Solar radiation pressure at Earth's distance from the sun (N/m^2)
    A = 1.0                   # Surface area of the spacecraft facing the sun (m^2)
    Cr = 1.3                  # Reflectivity coefficient
    sun_dir = normalize([1.0, 0.0, 0.0])  # Unit vector pointing from the sun

    # Compute the solar radiation pressure force:
    # F_srp = -P * A * Cr * sun direction
    F_srp = -P * A * Cr * sun_dir

    # Compute the torque from the solar force, assuming the force acts slightly off-center (0.1 m in x and y)
    T_srp = cross([0.1, 0.1, 0.0], F_srp)

    # Return both the SRP force and the resulting torque
    return F_srp, T_srp
end

# === Gravity Gradient Torque ===
# This function computes the torque on the spacecraft due to the gravity gradient effect.
# It happens because Earth's gravity pulls harder on one end of the spacecraft than the other.
#This function (Unless matrix is commented out) overides the inertial term in the original impax.yaml file
function compute_gravity_gradient(sim, r)
    u = sim.earth.mu   # Earth's gravitational parameter

    #+=+==+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    # ====!!!!! MANUAL INERTIA VALUE FOR EXPERIMENT!!!!! ====
    #Comment this out if user needs original inertial value from impax.yaml file

    #Cross product of original matrix caused gravity gradient to be zero so I
    #Attempted to play around with various values to view a potential gravity gradient
      I = [
        0.05  0.000001  0.000001;
        0.00000001  0.1  0.0000001;
        0.00000001  0.00000001  0.13
    ]

    R_hat = normalize(r)

    # Gravity gradient torque formula:
    # T = 3μ / |r|³ * (r̂ × (I * r̂))
    T_gg = (3 * u / norm(r)^3) * cross(R_hat, I * R_hat)

    return T_gg
end
# === Disturbance Type Classifier ===
#Perpahps issue with this if there's mixed forces... (I'm assuming there always will be?)
function classify_disturbance(torque)
    if norm(torque) < 1e-8
        return "Negligible"
    end
    dir = normalize(torque)

    if isapprox(dir, [1, 0, 0]; atol=0.2)
        return "Likely gravity gradient"
    elseif isapprox(dir, [0, 1, 0]; atol=0.2)
        return "Possible drag disturbance"
    elseif isapprox(dir, [0, 0, 1]; atol=0.2)
        return "Possible solar radiation pressure"
    else
        return "Unclassified or mixed"
    end
end

# === Main Disturbance Selector ===
# This function chooses the dominant disturbance type based on altitude and velocity.
function main_disturbance_type(r, v)
    altitude = norm(r) - 6371e3  # Approximate altitude by subtracting Earth's radius
    speed = norm(v)

    if altitude < 400e3
        return "Drag"
    elseif altitude < 800e3 && speed > 7500
        return "Solar"
    else
        return "Gravity Gradient"
    end
end

# === Main Disturbance ===
#This will display the disturbance forces and torque on the spacecraft
function compute_and_return_disturbances(sim, position_I, velocity_I, attitude_BI, time)
    println("=== Disturbance Force and Torque Analysis ===")

    torques = Vector{Vector{Float64}}()
    magnitudes = Float64[]

    for k in 1:length(time)
        r = position_I[:, k]
        v = velocity_I[:, k]
        dcm = attitude_BI[k]

        disturbance = main_disturbance_type(r, v)

        F_total = zeros(3)
        T_total = zeros(3)

        if disturbance == "Drag"
            F_total, T_total = compute_drag(r, v)
        elseif disturbance == "Solar"
            F_total, T_total = compute_solar_pressure()
        elseif disturbance == "Gravity Gradient"
            T_total = compute_gravity_gradient(sim, r)
        end

        push!(torques, T_total)
        push!(magnitudes, norm(T_total))

        # === Reduce Output Frequency (every 1000... This is an attempt to end the "endless loop")
    ### If we take average per 10000 time step of compute disturbances... maybe we can still get data
    ### Without the "endless loop"... Will the affect the accuracy of the data though?
        if k % 1000 == 1 || k == length(time)
            println("Time: $(round(time[k], digits=1)) s")
            println("  Total Torque: $(round.(T_total; digits=4)) Nm")
            println("  Likely Cause: $disturbance\n")
        end
    end

    # === Summary Output ===
    println("\n===== Disturbance Torque Summary =====")
    println("Mean torque magnitude: $(round(mean(magnitudes), digits=6)) Nm")
    println("Max torque magnitude:  $(round(maximum(magnitudes), digits=6)) Nm")
    println("Min torque magnitude:  $(round(minimum(magnitudes), digits=6)) Nm")

    return torques, magnitudes
end

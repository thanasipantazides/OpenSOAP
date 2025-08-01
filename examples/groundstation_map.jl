using OpenSOAP                     
using CairoMakie                   # For plotting and visualizing the passes
using Statistics                   # Min, max mean values
using FileIO                       # For loading image!
using SatelliteToolboxTransformations  # Earth transformations
import SatelliteToolboxTransformations: fetch_iers_eop  # Imports Earth orientation function (Eop)

# Function to find all continuous segments where the satellite is visible!
function pass_segments(vis::AbstractVector{<:Real})
    segs = Tuple{Int,Int}[]  # Stores start and end of visible segments
    in_pass = false          # Tracks whether we’re currently in a pass
    s = 0                    # Starts the index of current segment
    @inbounds for i in eachindex(vis)  # Iterates each visibility flag (1 or 0)
        if vis[i] == 1 && !in_pass     # Starts a new pass
            s = i
            in_pass = true
        elseif vis[i] == 0 && in_pass  # Ends the current pass
            push!(segs, (s, i - 1))    # Saves the pass segment
            in_pass = false
        end
    end
    if in_pass  # (Clumps the pass together until "pass zone" ends)
        push!(segs, (s, length(vis)))
    end
    return segs
end

# This function converts the satellite positions from ECI (inertial) to latitude and longitude
function positions_to_lla_deg(position_I, time, eops)
    lats = Vector{Float64}(undef, length(time))  # latitude array
    lons = Vector{Float64}(undef, length(time))  # longitude array
    @inbounds for k in eachindex(time)
        pos_lla = position_lla(position_I[:, k], time[k], eops)  # Convert!
        lats[k] = rad2deg(pos_lla[1])  # Convert from radians to degrees!
        lons[k] = rad2deg(pos_lla[2])
    end
    return lats, lons
end


#MAIN FUNCTION!
function plot_pass_shapes()
    println("Loading EOPs…")
    eops = fetch_iers_eop()  # Gets the Earth Orientation for transformations

    println("Loading mission…")
    sim = load_mission(joinpath("config", "mission.yaml"))  # Loads the mission config

    println("Integrating…")
    soln = integrate_system(dynamics_orbit!, sim.initstate, sim.tspan, sim.dt, sim)  # Runs the simulation!
    time = soln["time"]                   # This extracts the time vector
    Δt   = sim.dt                         # Time step
    position_I = soln["state"][1:3, :]    # Satellite position in ECI frame

    lats, lons = positions_to_lla_deg(position_I, time, eops)  # Converts to latitude/longitude

    ground_targets = [t for t in sim.mission.targets if t isa GroundTarget]  #Ground Stations :)

    println("Computing pass segments…")
    all_durations = Float64[]                       # This stores all pass durations
    per_station_segments = Dict{String, Vector{Tuple{Int,Int}}}()  # Segment map per station
    per_station_vis      = Dict{String, Vector{Float64}}()         # Visibility vector per station

    for gs in ground_targets
        vis = visibility_history(gs, soln)       # Computes the visibility vector for each station
        segs = pass_segments(vis)                # Segments the vector
        per_station_segments[gs.name] = segs     # Stores segments
        per_station_vis[gs.name] = vis           # Stores visibility flags

        for (s, e) in segs
            push!(all_durations, (e - s + 1) * Δt / 60)  # Stores the duration in minutes
        end
    end

    maxdur = isempty(all_durations) ? 1.0 : maximum(all_durations)  # Max duration for color scaling
    cmap = Makie.resample_cmap(:plasma, 256)  # This creates the color gradient (I chose a random # of colors(256))

    println("Plotting…")
    fig = Figure(size = (1400, 700), backgroundcolor = :white)  # Creates a new figure
    ax = Axis(fig[1, 1],
        title  = "Ground-Station Access Shapes",
        xlabel = "Longitude (°)",
        ylabel = "Latitude (°)",
        backgroundcolor = RGBf(1,1,1),
        xgridstyle = :dash, ygridstyle = :dash,
    )
    xlims!(ax, -180, 180)  # Sets the range for longitude
    ylims!(ax,  -90,  90)  # Sets the range for latitude

    # Optional: background world map
    try
        img = load("assets/map_diffuse.png")           # Loads image from assets file
        img = rotl90(img, 3)                            # This Rotates the image
        imgplot = image!(ax, -180..180, -90..90, img)   # Stretchs the image
        translate!(imgplot, 0, 0, -10)                   # Pushs image behind everything
    catch 
        @warn "Could not load world map image: $e"       # This prints if the image fails
    end
    lines!(ax, lons, lats, color = (:gray, 0.15), linewidth = 0.5)

    station_colors = Makie.wong_colors()  #Makes different colors for station mark

    for (j, gs) in enumerate(ground_targets)
        segs = per_station_segments[gs.name]  # Get segments for each station
        for (s, e) in segs
            dur_min = (e - s + 1) * Δt / 60.0              # The duration in minutes
            color_idx = clamp(round(Int, (dur_min / maxdur) * 255) + 1, 1, 256)  # Map to colormap
            c = cmap[color_idx]                           # This adds the color
            lw = 1 + 2 * (dur_min / maxdur)               # Line width increases with duration

            lines!(ax, lons[s:e], lats[s:e], color = c, linewidth = lw)  #This makes the "arcs" of the passes
        end

        lat_deg = rad2deg(gs.lla[1])  # Converts latitude to degrees
        lon_deg = rad2deg(gs.lla[2])  # Converts longitude to degrees

        scatter!(ax, [lon_deg], [lat_deg],
                 color = station_colors[(j-1) % length(station_colors) + 1],
                 markersize = 6)  # This shows the station location

        text!(ax, gs.name, position = (lon_deg + 2, lat_deg + 2), fontsize = 10, color = :black)  # Label for visual!
    end

    # This adds color bar to show duration meaning
    Colorbar(fig[1, 2], limits = (0, maxdur), colormap = cmap,
             label = "Pass duration (min)")

    display(fig)  #Shows the visual!!!
end

plot_pass_shapes()

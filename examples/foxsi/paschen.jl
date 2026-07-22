using Dates
using LinearAlgebra
using SatelliteToolboxAtmosphericModels
using CairoMakie

function paschen_vb(d, p; A=11.25, B=255.0, γ=0.02)
    return B*p*d / (log(A*p*d) - log(log(1+1/γ)))
end

function main(savepath::AbstractString)
    SpaceIndices.init()
    
    R_spec = 287.05
    h = 0:100:200e3
    lat = 65.1297*pi/180 # pfrr
    lon = -147.486*pi/180 # pfrr
    atm_all = AtmosphericModels.nrlmsise00.(DateTime("2024-04-17T22:13:00"), h, lat, lon)
    
    p = [atm.total_density*R_spec*atm.temperature for atm in atm_all]
    d = (0.1:0.005:10)*1e-3 # pitch grid
    # d = 10.0 .^(-4:1:-1)
    
    println("span: h: ", length(h), " d: ", length(d))
    P = [thisp for thisd in d, thisp in p]
    D = [thisd for thisd in d, thisp in p]
    H = [thish for thisd in d, thish in h]
    vb = [paschen_vb(thisd, thisp) for thisd in d, thisp in p]
    
    CairoMakie.activate!(type="pdf")
    fig = Figure()
    ax = Axis(
        fig[1,1],
        xlabel="Altitude [km]",
        ylabel="Electrode gap [mm]",
        title="Paschen-derived breakdown voltage profile\nPFRR, April 17, 2024",
        yscale=log10,
        yticks=([0.1,1.0,10.0], ["0.1", "1", "10"])
    )
    # xlims!(ax, [0,100])
    
    # levels = Makie.get_tickvalues(Makie.LinearTicks(7), extrema(vb)...)
    cpl = contourf!(ax, H / 1e3, D * 1e3, vb, levels=100:50:400, colormap=:rainbow)
    Colorbar(fig[1, 2], cpl, label="Breakdown voltage [V]")
    save(savepath, fig, transparent=true)
end
using OpenSOAP
import CSV
import Dates
import GLMakie
import CairoMakie
using GeometryBasics
using DataFrames

function main(outfile::AbstractString)

    MakieBackend = CairoMakie
    MakieBackend.activate!()
    fig = MakieBackend.Figure(size=(600,600))
    ax_batt = MakieBackend.Axis(fig[1,1], xlabel="Time [UTC]", ylabel="Battery [Wh]")
    ax_stor = MakieBackend.Axis(fig[2,1], xlabel="Time [UTC]", ylabel="Storage [MB]")

    infiles = [
        joinpath("log", "impax-inc55deg.csv"),
        joinpath("log", "impax-inc75deg.csv"),
        joinpath("log", "impax-inc95deg.csv")
    ]
    fnames = ["55 deg", "75 deg", "95 deg"]
    frames = []
    minlen = 1e12
    for f in infiles
        data = CSV.read(f, DataFrame; delim=',')
        if length(data[:,1]) < minlen
            minlen = length(data[:,1])
        end
        push!(frames, data)
    end
    for (k,f) in enumerate(frames)
        f = f[1:minlen, :]
        rename!(f, [:utc, :met, :t, :count, :pos_x, :pos_y, :pos_z, :vel_x, :vel_y, :vel_z, :att_11, :att_12, :att_13, :att_21, :att_22, :att_23, :att_31, :att_32, :att_33, :batt, :stor, :mode, :target])
        utcs = Dates.DateTime.(f[:, :utc], Dates.dateformat"dd u YYYY HH:MM:SS.sss")
        batts = f[:, :batt]
        MakieBackend.lines!(ax_batt, utcs, f[:, :batt] / 3600, label=fnames[k])
        MakieBackend.lines!(ax_stor, utcs, f[:, :stor] / 8e6, label=fnames[k])
    end
    
    fig[1:2,2] = MakieBackend.Legend(fig, ax_batt, "Inclination", framevisible=false)

    if MakieBackend === GLMakie
        display(fig)
    elseif MakieBackend === CairoMakie
        MakieBackend.save(outfile, fig, transparent=true)
    end

    # println(data[:,1])
    # targets = data[:,:target]
    # targets = data[:,:target]
    # # replace IDLE target with NaN for plotting
    # targets = [t == typemax(OpenSOAP.IDType) ? NaN : t for t in targets]

    
    
end
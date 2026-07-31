using Distributed
if nprocs() < 3
    addprocs(3-nprocs())
end

println("running with $(nprocs()) workers")

@everywhere import OpenSOAP

function main()
    fconfig = joinpath("config", "example.jsonc")
    println("launching core...")
    fcore = @spawnat workers()[1] OpenSOAP.run(config_path = fconfig)
    sleep(3)
    println("launching monitor...")
    # fmon = @spawnat workers()[2] OpenSOAP.show()
    OpenSOAP.monitor(config_path = fconfig)

    # f = fetch(fmon)
    f = fetch(fcore)

    # interrupt(fmon.where)
    sleep(2)
    interrupt(fcore.where)

    # # sleep(1)
    if isfile(OpenSOAP.unixsockname)
        rm(OpenSOAP.unixsockname)
    end
end

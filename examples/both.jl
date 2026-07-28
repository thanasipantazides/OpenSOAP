using Distributed
if nprocs() < 3
    addprocs(3-nprocs())
end

println("running with $(nprocs()) workers")

@everywhere import OpenSOAP

function main()

    println("launching core...")
    fcore = @spawnat workers()[1] OpenSOAP.run()
    sleep(3)
    println("launching monitor...")
    # fmon = @spawnat workers()[2] OpenSOAP.show()
    OpenSOAP.show()

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

using Distributed

if nprocs() < 3
    addprocs(3-nprocs())
end

function julia_main()::Cint

    fconfig = ARGS[1]
    println("launching core...")
    fcore = @spawnat workers()[1] run(config_path = fconfig)
    sleep(3)
    println("launching monitor...")
    fmon = @spawnat workers()[2] monitor(config_path = fconfig)
    # OpenSOAP.monitor(config_path = fconfig)

    fm = fetch(fmon)

    return 0
end

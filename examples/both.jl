using Distributed
import OpenSOAP



function main()

    if nprocs() < 3
        addprocs(3-nprocs())
    end
    
    println("running with $(nprocs()) workers")

    fmon = @async begin
        sleep(5) # wait for the core process to start
        println("launching monitor process")
        OpenSOAP.show()
    end
    
    println("launching core process")
    OpenSOAP.run()


    # note: with the two processes exchanging blocking like this, things work, but there is a low cap on how fast the sim can run. Around a time gain of 100:1.
end

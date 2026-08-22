
function launch()
    # eventually, calling this should:
    #   load up the config
    #   launch (in other processes) the core and monitor
    #   feed the core the config

    println("waiting for core to connect...")
    sock_repl = setup_server(SOAP_HOST, SOAP_REPL_PORT)
    println("connected!")
    return sock_repl
end

function poll()
    Threads.@spawn while true

    end
end

function update(sock::SocketWrapper, value::NetworkMessage)
    if value isa SimConfig
        error("Can't modify SimConfig while running!")
    end
    write_transport(sock, packetize(value, 0x0001))
end

using OpenSOAP

function make_solar_panels()::Vector{SolarPanel}
    pxf = SolarPanel([1.0;0;0], 0.295, 0.003018*6)
    pxb = SolarPanel([-1.0;0;0], 0.295, 0.003018*6)
    pzf = SolarPanel([0.0;0;1.0], 0.295, 0.003018*16)
    return [pxf; pxb; pzf]
end

function tumble_main()
    n = 10000
    s_I = [0;0;1]
    panels = make_solar_panels()
    A = zeros(n)
    for i in 1:n
        C_BI = r_random()
        for p in 1:length(panels)
            cosang = panels[p].normal'*C_BI*s_I
            if cosang > 0
                A[i] += panels[p].area*cosang
            end
        end
    end
    println(sum(A)/n)
end
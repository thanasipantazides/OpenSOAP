using Test
using OpenSOAP

@testset "Fixed-length data types" begin
    @testset "IDVector" begin

        # A normal vector of IDType
        ids = OpenSOAP.IDType[0x00, 0x01, 0x04, 0x478, typemax(OpenSOAP.IDType)]
        # Construct (by conversion) an IDVector, with explicit bucket size
        idv = OpenSOAP.IDVector{OpenSOAP.SOAP_MAX_ID_BUCKET}(ids)
        # Construct another IDVector with a different bucket size
        idp = OpenSOAP.IDVector{OpenSOAP.SOAP_MAX_ID_BUCKET ÷ 2}(ids)
        # Construct yet another, just to make sure a default bucket size is used
        idc = OpenSOAP.IDVector(ids)

        @test length(ids) == length(idv)
        @test all(ids .== idv)
        @test length(idv) == length(idp)
        @test all(idv .== idp)

        for k in eachindex(ids)
            @test ids[k] == idv[k]
            @test idv[k] == idp[k]
        end

    end

    @testset "FixedString" begin
        str1 = "Hello → 👍 unicode, λ, Malmö and Dragør!\n\t more"
        str2 = "He just kept talking in one long seemingly unbroken sentence so that no one had the chance to interrupt him it was really quite hypnotic. \n\nLet us go then, you and I, when the evening is spread out against the sky like a patient etherized upon the table. Let us go through half-deserted streets, the quick retreats..."

        @test_throws ArgumentError OpenSOAP.SizedString(str2)

        sstr1 = OpenSOAP.SizedString(str1)
        @test sstr1 == str1
        @test String(sstr1) == str1
    end

end

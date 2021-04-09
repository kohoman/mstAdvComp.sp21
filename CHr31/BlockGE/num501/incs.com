C ... common bge dimension parameter ...

        integer bnmax
        parameter(bnmax=512)

c ... common mg dimension parameter ...

        integer gptmax,npmax
        parameter(gptmax=360000,npmax=1000)       

C ... common bge dimension parameter ...

        integer bnmax
        parameter(bnmax=256)

c ... common mg dimension parameter ...

        integer gptmax,npmax
        parameter(gptmax=92500,npmax=260)       

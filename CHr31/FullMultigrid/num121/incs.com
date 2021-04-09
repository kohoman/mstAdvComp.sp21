c ... common dimension parameter ...

        integer gptmax,npmax
        parameter(gptmax=360000,npmax=1000)

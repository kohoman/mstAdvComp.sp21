c ... common dimension parameter ...

        integer gptmax,npmax
        parameter(gptmax=8000,npmax=200)

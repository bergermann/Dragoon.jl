



export unstuckDont, UnstuckDont



function unstuckDont(booster,hist,freqs,objFunction,x,f,args; showtrace=false)
    showtrace && println("No unstucking tried. Terminating.")

    return true
end

const UnstuckDont = Callback(unstuckDont)

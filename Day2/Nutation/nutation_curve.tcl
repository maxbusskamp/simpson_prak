spinsys {
    channels 1H
    nuclei 1H
    shift 1 0p 0p 0.0 0 0 0
}

par {
    # we can set an arbitrary RF field
    variable rf       50e3
    # and choose to increment the pulse duration by this amount in µs
    variable nutation_increment 1
    proton_frequency  300e6
    method            direct
    crystal_file      rep100
    spin_rate         10000
    gamma_angles      40
    sw                1e6/nutation_increment
    np                100
    start_operator    I1z
    detect_operator   I1p
    num_cores         64
    verbose           0100

}

proc pulseq {} {
    global par
    reset 
    # record the FID
    acq_block {
        # increment the time between acquisition points by a pulse
        # this cannot be done on a real spectrometer! why?
    	pulse $par(nutation_increment) $par(rf) y
    }
}

proc main {} {
    global par
    # make simulation
    set f [fsimpson]
    fsave $f nutation_rf_$par(rf).dat -xreim
    funload $f
}

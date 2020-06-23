#!/bin/bash

#Version 1
#This script builds a custom simpson input file for the respiration cp simulation. Custom Variables can be defined via the user input variables below.
#User Input Variable:
sim_name=$1
rf_power=$2
pulse_length=$3

rm $1.in
cat <<EOT >> $1.in
spinsys {
    channels 19F
    nuclei 19F
    shift 1 0p 10p 0.5 0 0 0 
}

par {
    method           direct
    proton_frequency 500e6
    spin_rate        100000
    crystal_file     rep20
    gamma_angles     10
    np               2048
    start_operator   I1z
    detect_operator  I1p
    sw               2000
    variable tsw     1e6/sw
    verbose          1101
	variable rf $rf_power
	variable p90 $pulse_length
}

proc pulseq {} {
    global par

    pulse \$par(p90) \$par(rf) y
    acq_block {
	delay \$par(tsw)
    }
}

proc main {} {
    global par
    
    set par(type) "real"

    set f [fsimpson]

    faddlb \$f 75 0

    # zerofilling to 2^N number
    fzerofill \$f 4096

    # make fft
    fft \$f

    # save spectrum in ACSII format
    fsave \$f \$par(name)_rf_\$par(p90).xy -xreim
    funload \$f
}




EOT

simpson $1.in
rm $1.in

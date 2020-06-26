spinsys {
    channels 13C
    nuclei 13C
    shift 1 0 0 0 0 0 0
}

par {
	proton_frequency	300e6
    spin_rate        	2500

	method				gcompute
    crystal_file     	rep256
    gamma_angles     	40
	
	sw               	spin_rate*gamma_angles
    variable tsw     	1e6/sw
    
	np               	4*1024
	variable lb			100
 
	start_operator   	I1x
    detect_operator  	I1p
}

proc pulseq {} {
    global par

    acq_block {
	    delay $par(tsw)
    }
}


proc main {} {
    global par
		
	# include the required packages
	# load packages
	lappend ::auto_path ~/simpson_prak/simpson/
	lappend ::auto_path ./opt 1.1
	if {![namespace exists opt]} {
        package require opt 1.1
        package require simpson-utils
        namespace import simpson-utils::*
    }
    
	
	# load experimental spectrum
	set par(exp) [fload 13C_MAS_exp.spe]
	
	
	
	# Define function for optimization
    opt::function rms

	
	# Define parameter for optimization according to name start_value, inc., lower_limit, upper_limit
	# names can be picked randomly but have to be assigned to the proper variables in proc rms (see below)!!!
	opt::newpar iso 0 100 -10 40
	opt::newpar csa 40 200 0 150
	opt::newpar eta 0.0 1 0 1

	# Scan parameter space
	opt::scan iso
	opt::scan csa
	opt::scan eta
	
	# Wrtie header for output
    puts " Progress		iso       csa      eta       rms"

	# start optimization of function "rms" with variables declared above
 	opt::minimize 1.0e-5
	
	# Save final spectrum
	rms 1

}

# Definition of rms 
proc rms {{save 0}} {
  global par

  # Set parameters using correct notation and in accordance with your choice of opt::newpar
  set f [fsimpson [list \
    [list shift_1_iso ${opt::iso}p]\
    [list shift_1_aniso ${opt::csa}p] \
    [list shift_1_eta $opt::eta] \
  ]]

  # add linebroadening, zerofilling and FT (Do not change!!!)
  faddlb $f $par(lb) 0
  fzerofill $f [expr 4*$par(np)]
  fft $f

  # normalize and scale correctly
  fautoscale $f $par(exp) -re
  
  # save test spectrum
  fsave $f $par(name)_test.xy -xreim
  
  # calculate RMS, compare and show results
  set rms [frms $f $par(exp) -re ]
  if {$save == 1} {
	  puts [format " \[%s\] %10.3f %10.3f %10.3f %10.3f" \
	     FINAL $opt::iso $opt::csa $opt::eta $rms]
    fsave $f $par(name)_final.xy -xreim
  } else {
	  puts -nonewline [format " \[%s\] %10.3f %10.3f %10.3f %10.3f\015" \
	     [progress] $opt::iso $opt::csa $opt::eta $rms]
  }
  flush stdout
  funload $f
  return $rms
}

proc progress {} {
  global par
  
  if ![info exists par(progress)] { set par(progress) -1 }
  incr par(progress)
  set str {"*---" \
           "-*--" \
           "--*-" \
           "---*" \
           "--*-" \
	   	     "-*--"}
  return [lindex $str [expr $par(progress)%6]]
}



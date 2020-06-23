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
	lappend ::auto_path /usr/local/tcl
	lappend ::auto_path ./opt 1.1
	if {![namespace exists opt]} {
        package require opt 1.1
        package require simpson-utils
        namespace import simpson-utils::*
    }
    
	
	# load experimental spectrum
	set par(exp) [fload 13C_MAS_exp.spe]
	
	
	# Define function for optimization


	
	# Define parameter for optimization according to name start_value, inc., lower_limit, upper_limit
	# names can be picked randomly but have to be assigned to the proper variables in proc rms (see below)!!!


	# Scan parameter space

	
	# Wrtie header for output
    puts " Progress		iso       csa      eta       rms"

	# start optimization of function "rms" with variables declared above
	
	# Save final spectrum

}

# Definition of rms 
proc rms {{save 0}} {
  global par

  # Set parameters using correct notation and in accordance with your choice of opt::newpar
  set f [fsimpson [list \

  ]]

  # add linebroadening, zerofilling and FT (Do not change!!!)
  faddlb $f $par(lb) 0
  fzerofill $f [expr 4*$par(np)]
  fft $f

  # normalize and scale correctly
  fautoscale $f $par(exp) -re
    
  # calculate RMS, compare and show results
  set rms [frms $f $par(exp) -re ]
  if {$save == 1} {
	  puts [format " \[%s\] %10.3f" \
	     FINAL ]
    fsave $f $par(name)_final.spe
  } else {
	  puts -nonewline [format " \[%s\] %10.3f\015" \
	     [progress] ]
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



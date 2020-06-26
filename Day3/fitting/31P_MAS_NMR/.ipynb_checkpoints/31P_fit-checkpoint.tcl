spinsys {
 
}

par {
    proton_frequency	300e6
    spin_rate        	10000
	variable noise  	0.005

	method				gcompute
    crystal_file     	rep256
    gamma_angles     	20
	
	sw               	spin_rate*gamma_angles
    variable tsw     	1e6/sw
    
	np               	4096
	variable lb			500
 
	start_operator   	Inx
    detect_operator  	Inp
 
	verbose          	0000
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
	set par(exp) [fload 31P_MAS_exp.spe]
	
	
	
	# Define function for optimization
    opt::function rms

	
	# Define parameter for optimization according to name start_value, inc., lower_limit, upper_limit
	# names can be picked randomly but have to be assigned to the proper variables in proc rms (see below)!!!
	
	# Wrtie header for output
    puts " Progress   ... "

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
	  puts [format " \[%s\] %10.3f\015" \
	     FINAL $rms]
    fsave $f $par(name)_final.spe
  } else {
	  puts -nonewline [format " \[%s\] %10.3f\015" \
	     [progress] $rms]
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
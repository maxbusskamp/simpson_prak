package require optimization
package require math::linearalgebra

package provide simpson-utils 1.1

namespace eval simpson-utils {
    namespace export \
	pardist \
	fselectregions \
	varian_pars \
	faddechoes \
	fixref \
	phasespectum \
	kval \
	funload-all \
	repeatfid \
	extfid  \
	ctemp \
	progressbar \
	progressbar2 \
	progressbar3 \
	progressbar4 \
	starttimer \
	stoptimer \
	savefit \
	saveopt \
	rmsscan \
	protonfreq \
	rmscalc \
	fbell \
	florentz\
	fpulse\
	fzeroim \
	fdcoffset \
	parametergrid
    
    global pi
    set pi [expr 4*atan(1)]

    proc parametergrid args {
	if {[llength $args] < 3} {
	    puts "Usage: parametergrid <start> <end> <stepsize>"
	    exit
	}
	set start [lindex $args 0]
	set end [lindex $args 1]
	set delta [lindex $args 2]
	if {[expr $end-$start] < 0} {
	    puts stderr "Error: parametergrid: <end> must be larger than <start>!"
	    exit
	}
	set npar [expr round(($end-$start)/$delta)]
	for {set i 0} {$i <= $npar } {incr i} {
	    lappend parlist [expr $start+$i*$delta]
	}
	return $parlist
    }
    
    proc fdcoffset args {
	if {[llength $args] < 1 || [llength $args] > 2} {
	    puts "Usage: fdcoffset <FID> ?<fraction, 0-1>?"
	    puts "       <fraction> is the part from the to be used for evaluating the DC offset."
	    puts "        defualt value is 0.5."
	    exit
	}
	set fid [lindex $args 0]
	
	if {[llength $args] == 1} {
	    set frac 0.5
	} else {
	    set frac [lindex $args 1]
	}
	
	set type [fget $fid -type]
	if {![string equal $type "fid"]} {
	    puts stderr "Error: fdcoffset may only be applied to FID data!" 
	    exit
	}
	if {[expr abs($frac)] > 1.0}  {
	    puts "Error: fdcoffset: <fraction> should be given between 0 and 1!"
	    exit
	}
	set np [fget $fid -np]
	set start [expr int($np-$frac*$np)+1] 
	
	# determine average value of the last $frac points of the FID
	# for real and imaginary part separatly
	set tmpre 0.0
	set tmpim 0.0
	for {set i $start} {$i <= $np} {incr i} {
	    set re [findex $fid $i -re]
	    set im [findex $fid $i -im]
	    set tmpre [expr $tmpre + $re]
	    set tmpim [expr $tmpim + $im]
	}
	set valre [expr $tmpre/($start)]
	set valim [expr $tmpim/($start)]
	# ... and substrac this value/offset from all point of the FID
	fexpr $fid [list \$re-$valre] [list \$im-$valim]
    }
    
    proc fzeroim args {
	set nargs 1
	if {[llength $args] < $nargs} {
	    puts "Usage: fzeroim <FID>"
	    exit
	}
	set fid [lindex $args 0]
	set type [fget $fid -type]
	if {![string equal $type "fid"]} {
	    puts stderr "Error: fzeroim may only be applied to FID data!" 
	    exit
	}
	fexpr $fid [list \$re] [list 0.0*\$im]
    }
    
    
    proc fbell args {
	set nargs 3
	if {[llength $args] < $nargs} {
	    puts "Usage: fbell <spectrum> <bell width, Hz> <bell center, Hz> ?options?"
	    puts "?options?  -save (saves a spectrum with the used Gaussian bell)"
	    puts "           -verb"
	    exit
	}
	set verb 0
	set save 0
	if {[llength $args] > $nargs} {
	    set options [lrange $args $nargs end]
	    foreach opt $options {
		switch -- $opt {
		    -save {
			set save 1
		    }
		    -verb {
			set verb 1
		    }		
		}
	    }
	}
	set f [lindex $args 0]
	set g [fdupzero $f]
	set w [lindex $args 1]
	set c [lindex $args 2]
	set np [fget $f -np]
	set sw [fget $f -sw]
	set dw [expr $sw/$np]
	if {$save} {
	    set gauss [fdupzero $f]
	    set max [lindex [fmaxheight $f] 0]
	}
	
	for {set i 1} {$i <= $np} {incr i} {
	    set x [expr -$sw/2.0+$dw*($i-1)]
	    set re [findex $f $i -re]
	    set im [findex $f $i -im]
	    set scale [gaussnorm $x $c $w]
	    if {$save} {
		fsetindex $gauss $i [expr $scale*$max] 0.0
	    }
	    fsetindex $g $i [expr $scale*$re] [expr $scale*$im]
	}
	if {$save} {
	    fautoscale $gauss $f
	    fsave $gauss gauss.spe
	    puts "fbell: Gauss-bell saved to gauss.spe..."
	    funload $gauss
	}
	funload $f
	return $g
    }
    
    
    proc gaussnorm {x m s} {
        set pi [expr 4*atan(1)]
        return [expr exp(-pow($x-$m,2)/(2.0*pow($s,2)))]
    }
    
proc fpulse args {
	set nargs 3
	if {[llength $args] < $nargs} {
	    puts "Usage: fsinc <spectrum> <pulse power, Hz> <flip angle, deg.> "
	    puts "applies excitation (intensity) profile for a given pulse"
	    puts "lit.: Conc. Mag. Res. DOI 10.1002/cmr.a"
		puts "tested only for pi/2 pulse!"
	    exit
	}
	set pi [expr 4*atan(1)]
	set verb 0
	set save 0
	set f [lindex $args 0]
	set g [fdupzero $f]
	set w [lindex $args 1]
	set angle [lindex $args 2]
	set radangle [expr $angle/360.0*2*$pi]
	set np [fget $f -np]
	set sw [fget $f -sw]
	set dw [expr $sw/$np]
	for {set i 1} {$i <= $np} {incr i} {
	    set x [expr -$sw/2.0+$dw*($i-1)]
	    set re [findex $f $i -re]
	    set im [findex $f $i -im]
	    set scale [pulsenorm $x $w $radangle]
	    if {$save} {
		fsetindex $sinc $i [expr $scale*$max] 0.0
	    }
	    fsetindex $g $i [expr $scale*$re] [expr $scale*$im]
	}
	funload $f
	return $g
    }
    
    
	proc pulsenorm {offset rf_power radangle} {
        set pi [expr 4*atan(1)]
		set offset1 [expr $offset*2*$pi]
		set rf_power1 [expr $rf_power*2*$pi]
		set rad1 [expr (abs($offset1)/$rf_power1)*0.0174532925199]
		set theta [expr atan(1/$rad1)]
		set alpha [expr $radangle*sqrt(1+(($offset1/$rf_power1)*($offset/$rf_power1)))]
		set radalpha [expr 0.0174532925199*$alpha]
		return [expr (sin($theta))*sqrt(((sin($alpha))*(sin($alpha))+(((cos($theta))*(1-(cos($theta)))))))]
     #   return [expr (sin($theta)*(180.0/$pi))*sqrt(((sin($alpha)*(180.0/$pi))*(sin($alpha)*(180.0/$pi))+(((cos($theta)*(180.0/$pi))*(1-(cos($theta)*(180.0/$pi)))))))]
    }
	
	proc florentz args {
	set nargs 3
	if {[llength $args] < $nargs} {
	    puts "Usage: fbell <spectrum> <bell width, Hz> <bell center> "
	    puts "applies a lorentz window function on freq. domain data "
		puts "see also 'fbell'"
	    exit
	}
	set verb 0
	set save 0
	if {[llength $args] > $nargs} {
	    set options [lrange $args $nargs end]
	    foreach opt $options {
		switch -- $opt {
		    -save {
			set save 1
		    }
		    -verb {
			set verb 1
		    }		
		}
	    }
	}
	set f [lindex $args 0]
	set g [fdupzero $f]
	set w [lindex $args 1]
	set c  [lindex $args 2]
	set np [fget $f -np]
	set sw [fget $f -sw]
	set dw [expr $sw/$np]
	for {set i 1} {$i <= $np} {incr i} {
	    set x [expr -$sw/2.0+$dw*($i-1)]
	    set re [findex $f $i -re]
	    set im [findex $f $i -im]
	    set scale [lorentznorm $x $c $w]
	   # if {$save} {
		#fsetindex $gauss $i [expr $mult*$scale*$max] 0.0
	    #}
	    fsetindex $g $i [expr $scale*$re] [expr $scale*$im]
	}
	funload $f
	return $g
    }
    
    
	proc lorentznorm {x m s} {
        set pi [expr 4*atan(1)]
        return [expr (1/$pi)*($s/($s*$s+($x-$m)*($x-$m)))]
	}
		
	
    proc protonfreq args {
	if {[llength $args] != 2} {
	    puts "Usage: protonfreq <Nuc> <resonance freq Hz>"
	    exit
	}
	set nuc1 [lindex $args 0]
	set res [lindex $args 1]
	set gamh [gamma 1H]
	set gamnuc1 [gamma $nuc1]
	set res [expr ($gamh/$gamnuc1)*$res]
	return $res
    }
    
    proc saveopt args {
	global par
	set verb 0
	#set curtime [clock format [clock seconds] -format "%Y-%m-%dT%H:%M:%S"]
	set curtime [exec date]
	
	if {[llength $args] < 1} {
	    puts "Usage: saveopt <file name> ?options?"
	    #puts "?options?: -verb"
	    puts "Saves final paramaters from optimization using the opt package"
	    puts "and SIMPSON simulation parameters to <file name>."
	    exit
	}
	
	if {[llength $args] >  1} {
	    set nopt 1
	    set options [lrange $args 1 end]
	    foreach opt $options {
		switch -- $opt {
		    -verb {
			set verb 1
		    }
		}
		incr nopt
	    }
	}
	set filename [lindex $args 0]
	set pars $opt::_internal::allparameters

	# save parameters to next file
	if {[file exists $filename]} {
	    set j 1
	    while {[file exists $filename.$j]} {
		incr j
	    }
	    puts "Warning: $filename exists! Old optimization results are copied to $filename.$j"
	    file copy -force $filename $filename.$j
	}
	set out [open $filename w]

	# print final parameters to file and check if it exists
	puts $out "Optimization using $par(name) finished on $curtime with the following parameters:"
	if {[info exists par(rms)]} {
	    puts $out "RMS: $par(rms)"
	}
	foreach p $pars {
	    puts $out [format "$p = %.5f" [set opt::$p]]
	}
	puts $out ""
	puts $out "Simulations were performed using the following SIMPSON parameters:"
	set simpars [lsort [array names par]]
	foreach s $simpars {
	    puts $out "$s  $par($s)"
	}
	close $out
    }
    
    proc fselectregions args {
	if {[llength $args]  < 3} {
	    puts "Usage: <desc> fselectregions <dataset> <spin rate, Hz> <peakwidth, Hz> ?options?"
	    puts "Returns a list with spectral regions usable for frms"
	    puts "?options?: -setup (use this option when setting up the regions of interest)"
	    puts "           -limits <left, Hz> <right ,Hz>"
	    puts "           -iso <val, Hz> (center of ssb pattern; defaults is 0.0)"
	    puts "           -remove <list with number(s)> (deletes regions with number(s), note that the list starts from 1)"
	    puts "           -odd  (odd number sidebands are returned)"
	    puts "           -even (even number sidebands are returned)" 
	    puts "           -ssborder <order> (selects <order> number of ssbs)"
	    puts "           -save (selected regions are saved to regions.spe)"
	    exit
	}
	set f [lindex $args 0]
	set ni [fget $f -ni]
	if {$ni > 1} {
	    puts "fselectregions: only works with 1D datasets!"
	    exit
	}
	set sw [fget $f -sw]
	set np [fget $f -np]
	set ref [fget $f -ref]
	set srate [expr double([lindex $args 1])]
	set pwidth [lindex $args 2]
	set ledge [expr $sw/2.0+$ref]
	set redge [expr $ledge-$sw]
	set sval $ledge
	set iso 0.0
	set ssbincr 1
	set odd 0
	set even 0
	set ssborder 0
	set limits 0
	set save 0
	set remove 0
	set setup 0
	if {[llength $args] > 3} {
	    set nopt 1
	    set options [lrange $args 3 end]
	    foreach opt $options {
		switch -- $opt {		    
		    -odd {
			set odd 1
			set ssbincr 2
		    }
		    -even {
			set even 1
			set ssbincr 2
		    }
		    -ssborder {
			set ssborder 1
			set ssbincr [lindex $options $nopt]
			if {![string is integer -strict $ssbincr]} {
			    puts "fselectregions: -ssborder, value must be integer!"
			    exit
			}
		    }
		    -remove {
			set remove 1
			set removelist [lsort -integer [lindex $options $nopt]]
		    }
		    -setup {
			set setup 1
			set save 1
		    }
		    -limits {
			set limits 1
			set ledge [lindex $options $nopt]
			set redge [lindex $options [expr $nopt+1]]			
			if {![string is double -strict $ledge] || ![string is double -strict $redge]} {
			    puts "fselectregions: -limits, given limits must be doubles!"
			    exit
			}
			if {$ledge < $redge} {
			    puts "Error: <left> must be larger than <right>!"
			    exit
			}
			set sval [expr round($ledge/$srate)*$srate]
		    }
		    -save {
			set save 1
		    }
		    -iso {
			set iso [lindex $options $nopt]
			if {![string is double -strict $iso]} {
			    puts "fselectregions: -iso, value must be double!"
			    exit
			}
		    }
		}
		incr nopt
	    }
	}
	set nssb [expr round(($ledge-$redge+$iso)/$srate)]
    	if {$odd} {
	    set nssb [expr $nssb-1]
	    set istart 1
	} else {
	    set istart 0
	}
	for {set i $nssb} {$i >= $istart} {incr i -$ssbincr} {
	    set center [expr $sval-$i*$srate+$iso]
	    set left [expr $center-$pwidth/2.0]
	    set right [expr $center+$pwidth/2.0]
	    if {$right < $redge || $left > $ledge} {
	    	continue
	    }
	    lappend plist [format "%.f %.f" $left $right]
	    lappend reglist [format "%.f %.f" $left $right]
	}
	if {$setup} {
	    puts "fselectregions: Number regions selected is [llength $reglist]"
	}
	if {$remove} {
	    # start of region list is 1! 
	    set rll [llength $reglist]
	    foreach e $removelist {
		set e [expr $rll-$e]
		if {$setup} {
		    puts "Removing region with index: $e"
		}
		set reglist [lreplace $reglist $e $e]
		set plist [lreplace $plist $e $e]
	    }
	}
	if {$save} {
	    set vs [expr double([lindex [fmaxheight $f] 0])/3.0]
	    set g [fcreate -type spe -np $np -sw $sw -ref $ref]
	    foreach reg $plist {
		set l [lindex $reg 0]
		set r [lindex $reg 1]
		set npend [expr int((($l-$ref)/$sw + 0.5)*$np) + 1]
		set npstart [expr int((($r-$ref)/$sw + 0.5)*$np) + 1]
		if {$npstart < 0 || $npstart > $np || $npend < 0 || $npend > $np} {
		    continue
		}
		for {set i $npend} {$i <= $npstart} {incr i} {
		    fsetindex $g $i $vs 0.0
		}
	    }
	    fsave $g regions.spe
	    if {$setup} {
		puts "fselectregions: selected regions are saved to regions.spe"
	    }
	    funload $g
	}
	if {$setup} {
	    puts "fselectregions: remove -setup when you have finished setting up regions"
	    exit
	}
	return $reglist
    }
    
    proc starttimer {} {
	global par
	set par(starttime) [clock seconds]
	set par(niter) 1
    }
    
    proc stoptimer {} {
	global par
	set par(stoptime) [clock seconds] 
	if {[info exists par(starttime)]} {
	    set elaptime [expr $par(stoptime)-$par(starttime)]
	    puts "Time for calculation was: $elaptime seconds"
	    #puts "Time for calculation was: [clock format $elaptime -format %M:%S]"
	    if {$par(niter) > 1} {
		puts "Number of iterations: $par(niter)"
	    }
	}
    }
    
    proc savefit args {
	global par
	if {[llength $args] < 1} {
	    puts "Usage: savefit <file name> ?options?"
	    puts "?options?: -verb"
	    exit
	}
	if {[llength $args] >  1} {
	    set nopt 1
	    set options [lrange $args 1 end]
	    foreach opt $options {
		switch -- $opt {
		    -verb {
			set verb 1
		    }
		}
		incr nopt
	    }
	}
	set pars $opt::_internal::allparameters
	set filename [lindex $args 0]
	if {$par(niter) < 2} {
	    puts "Each optimization step is saved to: $filename (use tail -f $filename to follow)"
	    set out [open $filename w]
	    puts -nonewline $out "RMS    "
	    foreach p $pars {
		puts -nonewline $out [format "%s    " $p]
	    }
	    puts $out ""
	} else {
	    set out [open $filename a]
	}
	# print parameters to file
	puts -nonewline $out [format "%.3f " $par(rms)]
	foreach p $pars {
	    puts -nonewline $out [format "%.3f " [set opt::$p]]
	}
	puts $out ""
	close $out
    }

    proc rmscalc args {
	global xrms yrms rmsstart rmsend
	set verb 0
	set range 0
	set sx 1.0
	set sy 1.0
	array unset xrms
	array unset yrms
	
	if {[llength $args] < 2} {
	    puts "Usage: rmscalc <parameter> <rms file> ?options?"
	    puts "?options?: -range <from> <to>"
	    puts "           -scalex <value>"
	    puts "           -scaley <value>"
	    puts "           -verb"
	    exit
	}
	if {[llength $args] >  2} {
	    set nopt 1
	    set options [lrange $args 2 end]
	    foreach opt $options {
		switch -- $opt {
		    -verb {
			set verb 1
		    }
		    -range {
			set range 1
			set from [lindex $options $nopt]
			set to [lindex $options [expr $nopt+1]]
			if {![string is double -strict $from] || ![string is double -strict $to]} {
			    puts "Error: -range <from> <to> (given: $from and $to) must be double!"
			    exit
			}
		    }
		    -scalex {
			set sx [lindex $options $nopt]
			if {![string is double -strict $sx]} {
			    puts "Error: -scalex <value> must be double!"
			    exit
			}
		    }
		    -scaley {
			set sy [lindex $options $nopt]
			if {![string is double -strict $sy]} {
			    puts "Error: -scaley <value> must be double!"
			    exit
			}
		    }
		}
		incr nopt
	    }
	}
	set p [lindex $args 0]
	set file [lindex $args 1]
	set fp [open $file r]
	set data [split [read $fp] \n]
	close $fp
	if {$range} {
	    set rmsstart [expr [lsearch [join $data] $from]/2]
	    set rmsend [expr [lsearch [join $data] $to]/2]
	    if {$rmsstart < 0 || $rmsend < 0} {
		puts "Error: -range <from> <to> ($from or $to) not found in rms data file!"
		exit
	    }
	} else {
	    set rmsstart 1
	    set rmsend [expr [llength $data]-2]
	}
	# read data...
	for {set i $rmsstart} {$i <= $rmsend} {incr i} {
	    set xrms($i) [expr [lindex $data $i 0]/$sx]
	    set yrms($i) [expr [lindex $data $i 1]/$sy]
	}
	set n [array size xrms]
	# calculate parameters...
	set x1 [calc1var x 1]
	set x2 [calc1var x 2]
	set x3 [calc1var x 3]
	set x4 [calc1var x 4]
	set y1 [calc1var y 1]
	set yx [calc2var y x 1 1]
	set yx2 [calc2var y x 1 2]

	set div [expr (pow($x2,3)+$n*$x3*$x3+$x1*$x1*$x4-$x2*(2*$x1*$x3+$n*$x4))]
	set a [expr ($x2*$x2*$y1-$x1*$x3*$y1+$n*$x3*$yx+$x1*$x1*$yx2-$x2*($x1*$yx+$n*$yx2))/$div]
	set b [expr ($x1*$x4*$y1+$x2*$x2*$yx-$n*$x4*$yx+$n*$x3*$yx2-$x2*($x3*$y1+$x1*$yx2))/$div]
	set c [expr ($x3*$x3*$y1-$x2*$x4*$y1+$x1*$x4*$yx+$x2*$x2*$yx2-$x3*($x2*$yx+$x1*$yx2))/$div]
	if {$a < 0} {
	    puts "Error: a is below zero! Check that rms data includes a minimum (in given range?)"
	    exit
	}
	set min [expr -$b/(2*$a)]
	set dpar [expr 2.0/sqrt($a)]

	if {$verb} {
	    puts [format "$p: %7.2f +/- %5.4f" $min $dpar]
	    puts "Polynomium: f(x)=a*x**2+b*x+c"	    
	    puts "a=$a;b=$b;c=$c"
	    puts "$n points used in calculation"
	}
	return [list $min $dpar $a $b $c]
    }

    proc calc1var {par order} {
	global xrms yrms rmsstart rmsend
	set j 0
	if {[string equal $par x]} {
	    for {set i $rmsstart} {$i <= $rmsend} {incr i} {
		set j [expr $j+pow($xrms($i),$order)]
	    }
	} elseif {[string equal $par y]} {
	    for {set i $rmsstart} {$i <= $rmsend} {incr i} {
		set j [expr $j+pow($yrms($i),$order)]
	    }
	}
	return $j
    }
    
    proc calc2var {par1 par2 order1 order2} {
	global xrms yrms rmsstart rmsend
	set j 0
	if {[string equal $par1 x] || [string equal $par2 y]} {
	    for {set i $rmsstart} {$i <= $rmsend} {incr i} {
		set j [expr $j+pow($xrms($i),$order1)*pow($yrms($i),$order2)]
	    }
	} elseif {[string equal $par1 y] || [string equal $par2 x]} {
	    for {set i $rmsstart} {$i <= $rmsend} {incr i} {
		set j [expr $j+pow($yrms($i),$order1)*pow($xrms($i),$order2)]
	    }
	}
	return $j
    }
    
    proc rmsscan args {
	global par
	if {[llength $args] < 4} {
	    puts "Usage: rmsscan <parameter> <lower limit> <upper limit> <number>"
		puts "		 <parameter> is without opt::"
		puts "		 <number> is number of steps in each direction!"
	    exit
	}
	
	set par(verb) 0
	set par(progress) 0
	set verb 0
	
	#set parfile 0
	
	set parameter [lindex $args 0]
	set high [lindex $args 2]
	set low [lindex $args 1]
	if {$high <= $low} {
	puts "upper limit is below lower limit!"
	exit
	}
	set number [lindex $args 3]
	set delta [expr ($high-$low)/($number-1)]
	set parsave [set opt::${parameter}]
	if [info exists rms_data_list] {unset rms_data_list}
		
		
	#start parameter sweep
	set opt::${parameter} $low
		
	for {set i 0} {$i < $number} {incr i} {
		set p [set opt::${parameter}]
		set y [rms]
		set opt::${parameter} [expr $low+$i*$delta]
		lappend rms_data_list [list $p $y]
	}
	set opt::${parameter} $parsave
		
	set length [llength $rms_data_list]
	set sx4  0.0
	set sx3  0.0
	set sx2  0.0
	set sx   0.0
	set sx2y 0.0
	set sxy  0.0
	set sy   0.0
	for {set i 0} {$i < $length} {incr i} {
		set p [lindex $rms_data_list $i 0]
		set y [lindex $rms_data_list $i 1]
		set sx   [expr $sx + $p]
		set sx2  [expr $sx2 + $p*$p]
		set sx3  [expr $sx3 + $p*$p*$p]
		set sx4  [expr $sx4 + $p*$p*$p*$p]
		set sy   [expr $sy + $y]
		set sxy  [expr $sxy + $p*$y]
		set sx2y [expr $sx2y + $p*$p*$y]
	}

	set mat [::math::linearalgebra::mkMatrix 3 3 0]
	set vec [::math::linearalgebra::mkVector 3 0]
	::math::linearalgebra::setelem mat 0 0 $sx4
	::math::linearalgebra::setelem mat 0 1 $sx3
	::math::linearalgebra::setelem mat 0 2 $sx2
	::math::linearalgebra::setelem mat 1 0 $sx3
	::math::linearalgebra::setelem mat 1 1 $sx2
	::math::linearalgebra::setelem mat 1 2 $sx
	::math::linearalgebra::setelem mat 2 0 $sx2
	::math::linearalgebra::setelem mat 2 1 $sx
	::math::linearalgebra::setelem mat 2 2 $length

	::math::linearalgebra::setelem vec 0 $sx2y
	::math::linearalgebra::setelem vec 1 $sxy
	::math::linearalgebra::setelem vec 2 $sy

	set res [::math::linearalgebra::solveGauss $mat $vec]

	set a [lindex $res 0]
	if {$a > 0} {
		set answer [expr 2.0/sqrt($a)]
	}
	
	puts -nonewline "\n"
		    puts "RMS scan finished..."
			puts [format "%10.3f" $answer]
    }

	
	proc progress {} {
	global par
  
	if ![info exists par(progress)] { set par(progress) -1 }
	incr par(progress)
	set str {"*   " \
           " *  " \
           "  * " \
           "   *" \
           "  * " \
		   " *  "}
	return [lindex $str [expr $par(progress)%6]]
}
	
	
	
	
    proc progressbar {} {
	global par
	incr par(progress)
	switch [expr $par(progress)%8] {
	    0 {return "\[.    \]"}
	    1 -
	    7 {return "\[ o   \]"}
	    2 -
	    6 {return "\[  0  \]"}
	    3 -
	    5 {return "\[   o \]"}
	    4 {return "\[    .\]"}
	}
    }

    proc progressbar4 {} {
	global par
	incr par(progress)
	switch [expr $par(progress)%8] {
	    0 {return "."}
	    1 -
	    7 {return "o"}
	    2 -
	    6 {return "0"}
	    3 -
	    5 {return "o"}
	    4 {return "."}
	}
    }
    
    proc progressbar3 {} {
	global par
	incr par(progress)
	switch [expr $par(progress)%6] {
	    0 {return "-"}
	    1 -
	    6 {return "\\"}
	    2 -
	    3 {return "|"}
	    4 -
	    5 {return "/"}
	}
    }

    proc progressbar2 {} {
	global par
	incr par(progress)
	switch [expr $par(progress) % 18] {
	    0 {return "\[Z         \]"}
	    1 -
	    18 {return "\[ z        \]"}
	    2 -
	    17 {return "\[  Z       \]"}
	    3 -
	    16 {return "\[   z      \]"}
	    4 -
	    15 {return "\[    Z     \]"}
	    5 -
	    14 {return "\[     z    \]"}
	    6 -
	    13 {return "\[      Z   \]"}
	    7 -
	    12 {return "\[       z  \]"}
	    8 -
	    11 {return "\[        Z \]"}
	    9 -
	    10 {return "\[         z\]"}
	}
    }

    proc extfid {f ni_start ni_end mult lb1 gbfrac nz {zero_im 1}} {
	global pi
	set sw [fget $f -sw]
	set np [fget $f -np]
	set sw1 [fget $f -sw1]
	set dw1 [expr 1.0/$sw1]
	set ref [fget $f -ref]    
	set ref1 0.0
	set ni [expr $ni_end-$ni_start+1]
	set newni [expr $mult*$ni]
	set g [fcreate -type spe -np $np -ni $nz -sw $sw -sw1 $sw1 -ref $ref -ref1 $ref1]
	fzero $g
	puts "t1-dimension is extended from: $ni to $newni points"
	for {set  i 1} {$i <= $np} {incr i} {
	    for {set j $ni_start} {$j <= $ni_end} {incr j} {
		set re [findex $f $i $j -re]
		if {$zero_im} {
		set im 0.0
		} else {
		    set im [findex $f $i $j -im]
		}
		for {set k 0} {$k < $mult} {incr k} {
		    fsetindex $g $i [expr ($j-$ni_start+1)+$k*$ni] $re $im
		}
	    }
	}
	# line broadening...
	set lb [expr $lb1*(1.0-$gbfrac)]
	set gb [expr $lb1*$gbfrac]
	set fac [expr -$pi*$lb1*$dw1]
	set fag [expr 2.0*$pi*$gb*$dw1/(4.0*sqrt(log(2.0)))]
	puts "Applying linebroadening..."
	for {set i 1} {$i <= $np} {incr i} {
	    for {set j 1 } {$j <= $newni} {incr j} {
		set val [findex $g $i $j -re]
		set re [expr $val*exp($fac*$j-pow($fag,2))]
		set im [findex $g $i $j -im]
		fsetindex $g $i $j $re $im
	    }
	}
	return $g
    }
        
    proc repeatfid {f mult {zero_im 0}} {
	set np [fget $f -np]
	set sw [fget $f -sw]
	set newnp [expr $mult*$np]
	set g [fcreate -np $newnp -sw $sw]
	for {set i 1} {$i <= $np} {incr i} {
	    set re [findex $f $i -re]
	    if {$zero_im} {
		set im 0.0
	    } else {
		set im [findex $f $i -im]
	    }
        for {set j 1} {$j <= $mult} {incr j} {
            fsetindex $g [expr ($i+($j-1)*$np)] $re $im
        }
    }
	return $g
    }
    
    proc pardist args {
	if {[llength $args] < 2} {
	    puts "Usage: pardist <model> ?options?"
	    puts "<model>: normal      <mean> <sigma> <\#points>"
	    puts "       : half-normal <mean> <sigma> <\#points>"
	    puts "       : log-normal  <mean> <sigma> <\#points>"
	    puts "?options?: -plot (points are saved to pardist.plot and pardist.function, viewable with gnuplot)"
	    puts "           -verb"
	    exit
	}
	set verb 0
	set upto 0
	set limits 0
	set plot 0
	set model [lindex $args 0]
	switch -- $model {
	    normal {
		if {[llength $args] < 4} {
		    puts "pardist: $model, not enough input parameters given!"
		    puts "normal <mean> <sigma> <\#points>"
		    exit
		}
		set mean [lindex $args 1]
		set sigma [lindex $args 2]
		if {$sigma < 0.0} {
		    puts "pardist: $model, <sigma> must be larger than 0."
		    exit
		}
		set npoints [lindex $args 3]
		set ninput 4
	    }	    
	    half-normal {
		if {[llength $args] < 4} {
		    puts "pardist: $model, not enough input parameters given!"
		    puts "$model <value> <sigma> <\#points>"
		    exit
		}
		set mean [lindex $args 1]
		set sigma [lindex $args 2]
		set npoints [lindex $args 3]
		set ninput 4
	    }
	    log-normal {
		if {[llength $args] < 4} {
		    puts "pardist: $model, not enough input parameters given!"
		    puts "$model <mean value> <sigma> <\#points>"
		    exit
		}
		set mean [lindex $args 1]
		set sigma [lindex $args 2]
		set npoints [lindex $args 3]
		set ninput 4
		# if {$sigma < 0.0} {
		#puts "Error: pardist normal, <sigma> must be larger than 0."
		#exit
		#}
	    }
	    default {
		puts "pardist: $model is not a valid model!"
		puts "<model>: normal     <mean value> <sigma> <\#points>"
		puts "       : log-normal <mean value> <sigma> <\#points>"
		exit
	    }
	}
	if {$npoints < 2} {
	    return [list [list [expr double($mean)] 1.0]]
	}
	if {[llength $args] > $ninput} {
	    set nopt 1
	    set options [lrange $args $ninput end]
	    foreach opt $options {
		switch -- $opt {		    
		    -plot {
			set plot 1
		    }
		    -verb {
			set verb 1
		    }
		}
		incr nopt
	    }
	}
	
	# calculate weight depending on model
	switch -- $model {
	    normal {
		# check if number of points is odd
		set function "f(x)=1.0/(s*sqrt(2.0*pi))*exp(-(x-m)**2/(2*s**2))\nm=$mean;s=$sigma"
		if {![expr $npoints % 2]} {
		    puts "pardist: <\#points> must be odd with <model> = $model!" 
		    exit
		}		
		# make list with values
		set min [expr $mean-3.0*$sigma]
		set delta [expr ($mean-$min)/(($npoints-1)/2)]
		set delta2 [expr $delta/2.0]
		for {set i 0} {$i < $npoints} {incr i} {
		    set mval [expr $min+$i*$delta]
		    lappend lval [list $mval [expr $mval-$delta2] [expr $mval+$delta2]]
		}
		foreach v $lval {
		    set w [weight gauss $mean $sigma [lindex $v 1] [lindex $v 2]]
		    lappend res [list [format "%.2f" [lindex $v 0]] $w]
		}
	    }
	    half-normal {
		set sigma2 [expr double($sigma)/double($mean)]
		#set function "f(x)=sqrt(2.0/(pi*s**2))*exp(-x**2/(2*s**2))\ns=$sigma2"
		set function "f(x)=exp(-x**2/(2*s**2))\ns=$sigma2"
		set min [expr 3.0*$sigma2]
		set delta [expr $min/$npoints]
		set delta2 [expr $delta/2.0]
		for {set i 0} {$i < $npoints} {incr i} {
		    set mval [expr $i*$delta]
		    set w [halfgauss $mval $mean $sigma2]
		    lappend res [list [format "%.2f" [expr $mean-$mval*$mean]] [format "%6.4f" $w]]
		}
	    }
	    log-normal {
		puts "Not finished..."
		exit 
		set w [weight loggauss $mean $sigma [lindex $v 1] [lindex $v 2]]
		set function "f(x)=1.0/(x*s*sqrt(2.0*pi))*exp(-(log(x)-m)**2/(2*s**2))\nm=$mean;s=$sigma"
	    }
	}
	
	# make plot files for gnuplot
	if {$plot} {
	    set plotout [open pardist.function w]
	    set plotpoints [open pardist.points w]
	    puts $plotout $function
	    if {$model != "half-normal"} {
		puts $plotout "set xrange \[[expr $mean-5*$sigma]:[expr $mean+5*$sigma]\];"
	    } else {
		puts $plotout "set xrange \[0:1\];"
	    }
	    puts $plotpoints "\# data points for $model with mean = $mean and sigma = $sigma"
	    foreach p $res {
		puts $plotpoints "[lindex $p 0] [lindex $p 1]"
	    }
	    puts "$model distribution points saved to pardist.points and pardist.function"
	    close $plotpoints
	    close $plotout
	}
	return $res
    }
    
    proc weight {func mean std from to {nmax 100}} {
	# Integration is performed using Simpson's rule and split in to n composite intervals (default is 100).
	# See forexample: http://en.wikipedia.org/wiki/Simpson's_rule
	set y 0.0
	set delta [expr ($to-$from)/$nmax]
	for {set n 0} {$n < $nmax} {incr n} {
	    set x1 [expr $from+$n*$delta]
	    set x2 [expr $from+($n+1)*$delta]    
	    set y1 [$func $x1 $mean $std]
	    set y2 [$func [expr ($x1+$x2)/2.0] $mean $std]
	    set y3 [$func $x2 $mean $std]
	    set y [expr $y+($x2-$x1)/6.0*($y1+4*$y2+$y3)]
	}
	return [format "%6.4f" $y]
    }
    
    proc loggauss {x m s} {
	set pi [expr 4*atan(1)]
	return [expr 1.0/($x*$s*sqrt(2.0*$pi))*exp(-(log($x)-pow($m,2))/(2.0*pow($s,2)))]
    }

    proc gauss {x m s} {
	set pi [expr 4*atan(1)]
	return [expr 1.0/($s*sqrt(2.0*$pi))*exp(-pow($x-$m,2)/(2.0*pow($s,2)))]
    }
    
    proc halfgauss {x m s} {
	set pi [expr 4*atan(1)]
	#return [expr sqrt(2.0/($pi*pow($s,2)))*exp(-pow($x,2)/(2.0*pow($s,2)))]
	return [expr exp(-pow($x,2)/(2.0*pow($s,2)))]
    }
    

    proc funload-all {spec} {
	foreach s $spec {
	    funload $s
	}
    }
    
    proc varian_pars {base pars} {
	global specpar
	set par $base/procpar
	set fp [open $par r]
	set tmp [split [read $fp] \n]
	foreach p $pars {
	    set i [lsearch -exact -regexp $tmp "^$p "]	
	    set pv [lindex $tmp [expr $i+1] 1]
	    set specpar(${p}) $pv
	}
    }
 
    proc faddechoes {f g} {
	set fer [fget $f -ref]
	set fer1 [fget $f -ref1]
	fset $f -ref 0
	fset $g -ref 0
	fset $f -ref1 0
	fset $g -ref1 0
	fadd $f $g
	fset $f -ref $fer
	fset $f -ref1 $fer1
    }

    proc phasespectum {f lb nz rp lp} {
	fset $f -ni 0
	fset $f -ref 0
	faddlb $f $lb 0
	fzerofill $f $nz
	fft $f
	fphase $f -rp $rp -lp $lp
	fsave $f spec-1d.spe
	funload $f
	exec simplot spec-1d.spe &
	exit
    }
    
    proc fixref {file spectrum} {
	global specpar
	puts "exec sethdr $file -xCAR [expr [fget $spectrum -ref]/$specpar(sfrq)] -yCAR [expr ([fget $spectrum -ref1])/$specpar(sfrq)] -xOBS $specpar(sfrq) -yOBS $specpar(sfrq)"
	exec sethdr $file -xCAR [expr [fget $spectrum -ref]/$specpar(sfrq)] -yCAR [expr ([fget $spectrum -ref1])/$specpar(sfrq)] -xOBS $specpar(sfrq) -yOBS $specpar(sfrq)
    }
    
    proc kval {ival coherence} {
	switch -exact $ival {
	    3/2 {
		switch -exact $coherence {
		    3Q {
			set k [expr 7.0/9.0]
		    }
		    default {
			puts stderr "kval: Only 3Q possible for I = 3/2"
			exit
		    }
		}
	    }
	    5/2 {
		switch -exact $coherence {
		    3Q {
			set k [expr 19.0/12.0]
		    }
		    5Q {
			set k [expr 25.0/12.0]
		    }
		    default {
			puts stderr "kval: Only 3Q and 5Q"
			exit
		    }
		}
	    }	 
	    7/2 {
		switch -exact $coherence {
		    3Q {
			set k [expr 101.0/45.0]
		    }
		    5Q {
			set k [expr 11.0/9.0]
		    }
		    default {
			puts stderr "kval: Only 3Q and 5Q"
			exit
		    }
		}
	    }
	    9/2 {
		switch -exact $coherence {
		    3Q {
			set k [expr 91.0/36.0]
		    }
		    5Q {
			set k [expr 95.0/36.0]
		    }
		    default {
			puts stderr "kval: Only 3Q and 5Q"
			exit
		    }
		}   
	    }
	    default {
		puts stderr "kval:\nIval must be: 3/2, 5/2, 7/2, or 9/2"
		puts stderr "Coherence: 3Q and 5Q only"
		exit
	    }
	}
	return $k
    }

    proc ctemp {settemp vr} {
	#numbers in the equation have been determined by Ingo (a while ago...)
	#T(actual)=-12.1+1.019*T(displayed)+0.62*(spinning)-0.0023*T(displayed)*
	#spinning+0.0364*spinning^2
	# vr in kHz and T actual in K.	
	return [expr (-12.1+1.019*$settemp+0.62*$vr-0.0023*$settemp*$vr+0.0364*pow($vr,2))-273.15]
    }
}

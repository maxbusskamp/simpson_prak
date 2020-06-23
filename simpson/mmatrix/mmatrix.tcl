package provide mmatrix 1.1

namespace eval mmatrix {
    namespace export \
	mputmatrix \
	mgetmatrix \
	mdirprod \
	mmult \
	mdim \
	madd \
	msub \
	madjoint \
 	mtrace \
	mscale \
	mtranspose \
	mconjugate \
	mcreatematrix
    
    proc mmult {m1 m2} {
	set n [llength $m1]
	set m $m1
	for {set r 0} {$r < $n} {incr r} {
	    for {set c 0} {$c < $n} {incr c} {
		lset m $r $c {0 0}
		for {set i 0} {$i < $n} {incr i} {
		    lset m $r $c 0 [expr [lindex $m $r $c 0] \
					+ [lindex $m1 $r $i 0]*[lindex $m2 $i $c 0] \
					- [lindex $m1 $r $i 1]*[lindex $m2 $i $c 1]]
		    lset m $r $c 1 [expr [lindex $m $r $c 1] \
					+ [lindex $m1 $r $i 0]*[lindex $m2 $i $c 1] \
					+ [lindex $m1 $r $i 1]*[lindex $m2 $i $c 0]]
		}
	    }
	}
	return $m
    }

    proc madd {m1 m2} {
	if {[mdim $m1] != [mdim $m2]} {
	    puts stderr "Error: madd: Matrices must have same dimension"
	    exit
	}
	set n [llength $m1]
	set m $m1	
	for {set r 0} {$r < $n} {incr r} {
	    for {set c 0} {$c < $n} {incr c} {
		lset m $r $c 0 [expr [lindex $m1 $r $c 0]+[lindex $m2 $r $c 0]]
		lset m $r $c 1 [expr [lindex $m1 $r $c 1]+[lindex $m2 $r $c 1]]
	    }
	}
	return $m
    }


    proc msub {m1 m2} {
	if {[mdim $m1] != [mdim $m2]} {
	    puts stderr "Error: msub: Matrices must have same dimension"
	    exit
	}
	set n [llength $m1]
	set m $m1
	for {set r 0} {$r < $n} {incr r} {
	    for {set c 0} {$c < $n} {incr c} {
		lset m $r $c 0 [expr [lindex $m1 $r $c 0]-[lindex $m2 $r $c 0]]
		lset m $r $c 1 [expr [lindex $m1 $r $c 1]-[lindex $m2 $r $c 1]]
	    }
	}
	return $m
    }
    
    proc madjoint {m} {
	set n [llength $m]
	set m1 $m
	for {set r 0} {$r < $n} {incr r} {
	    for {set c 0} {$c < $n} {incr c} {
		lset m1 $r $c [list [lindex $m $c $r 0] [expr -[lindex $m $c $r 1]]]
	    }
	}
	return $m1
    }

    proc mtranspose {m} {
	set n [llength $m]
	set m1 $m
	for {set r 0} {$r < $n} {incr r} {
	    for {set c 0} {$c < $n} {incr c} {
		lset m1 $r $c [list [lindex $m $c $r 0] [expr [lindex $m $c $r 1]]]
	    }
	}
	return $m1
    }

    proc mconjugate {m} {
	set n [llength $m]
	set m1 $m
	for {set r 0} {$r < $n} {incr r} {
	    for {set c 0} {$c < $n} {incr c} {
		lset m1 $r $c 1 [expr -[lindex $m $r $c 1]]
	    }
	}
	return $m1
    }
    
    proc mscale {m val} {
	set n [llength $m]
	set m1 $m
	for {set r 0} {$r < $n} {incr r} {
	    for {set c 0} {$c < $n} {incr c} {
		lset m1 $r $c 0 [expr [lindex $m $r $c 0]*$val]
		lset m1 $r $c 1 [expr [lindex $m $r $c 1]*$val]
	    }
	}
	return $m1
    }

    proc mtrace {m} {
	set re 0.0
	set im 0.0
	set l [llength $m]
	for {set i 0} {$i < $l} {incr i} {
	    set re [expr $re + [lindex $m $i $i 0]]
	    set im [expr $im + [lindex $m $i $i 1]]
	}
	return [list $re $im]
    }
    
    proc mputmatrix {m {part "reim"} {fm "%4.2g"}} {
	foreach i $m {
	    foreach j $i {
		if {[llength $j] == 2} {
		    if {$part == "reim"} {
			puts -nonewline [format "($fm,$fm) " [lindex $j 0] [lindex $j 1]]
		    }
		    if {$part == "re"} {
			puts -nonewline [format "$fm " [lindex $j 0]]
		    }
		    if {$part == "im"} {
			puts -nonewline [format "$fm " [lindex $j 1]]
		    }
		} else {
		    puts -nonewline [format $fm $j]
		}
	    }
	    puts ""
	}
    }
    
    proc mcreatematrix {r c} {
	for {set i 0} {$i < $r} {incr i} {
	    append m "{"
	    for {set j 0} {$j < $c} {incr j} {
		append m "{0 0} "
	    }
	    set m [string trimright $m]
	    append m "} "
	}
	return $m
    }
    
    proc mdim {m} {
	return [llength [lindex $m 0]]
    }

    proc mgetmatrix {spin op} {
	global par
	for {set i 2} {$i <= 24} {incr i} {
	    set sq${i} [expr sqrt($i)]
	    set sq${i}h [expr sqrt($i)/2]
	}

	# create empty matrix based on spin number
	# lset m row col real/im val
	switch -exact -- $spin {
	    1/2 {set m [mcreatematrix 2 2]
		set s 0.5}
	    1   {set m [mcreatematrix 3 3]
		set s 1.0}
	    3/2 {set m [mcreatematrix 4 4] 
		set s 1.5}
	    5/2 {set m [mcreatematrix 6 6]
		set s 2.5}
	    7/2 {set m [mcreatematrix 8 8]
		set s 3.5}
	    9/2 {set m [mcreatematrix 10 10]
		set s 4.5}
	    default {
		puts "Usage: mgetmatrix <spin> <operator>"
		puts "Wrong spin number"
		puts "spin should be: 1/2, 1, 3/2, 5/2, 7/2, or 9/2"
		exit
	    }
	}

	# change appropriate elements based on spin number
	switch -exact -- $spin {
	    1/2 {
		switch -exact -- $op {
		    E {
			lset m 0 0 0 1
			lset m 1 1 0 1
		    }
		    Ix {
			lset m 0 1 {0.5 0}
			lset m 1 0 {0.5 0}
		    }
		    Iy {
			lset m 0 1 {0 -0.5}
			lset m 1 0 {0 0.5}
		    }
		    Iz {
			lset m 0 0 {0.5 0}
			lset m 1 1 {-0.5 0}
		    }
		    Ip {
			lset m 0 1 {0.5 0}
		    }
		    Im {
			lset m 1 0 {0.5 0}
		    }
		    default {
			puts "Usage: mgetmatrix <spin> <operator>"
			puts "Not an operator"
			puts "Options are: E Ix Iy Iz Ip Im"
			exit
		    }
		}
	    }
	    1 {
		switch -exact -- $op {
		    E {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 1
			}
		    }
		    Iz {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 [expr $s-$i]
			}
		    }
		    Ix {
			lset m 0 1 0 $sq2h
			lset m 1 0 0 $sq2h
			lset m 1 2 0 $sq2h
			lset m 2 1 0 $sq2h
		    }
		    Iy {
			lset m 0 1 1 -$sq2h
			lset m 1 0 1 $sq2h
			lset m 1 2 1 -$sq2h
			lset m 2 1 1 $sq2h
		    }
		    Ip {
			lset m 0 1 0 $sq2
			lset m 1 2 0 $sq2
		    }
		    Im {
			lset m 1 0 0 $sq2
			lset m 2 1 0 $sq2
		    }
		}
	    }
	    3/2 {
		switch -exact -- $op {
		    E {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 1
			}
		    }
		    Iz {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 [expr $s-$i]
			}
		    }
		    Ix {
			lset m 0 1 0 $sq3h
			lset m 1 0 0 $sq3h
			lset m 1 2 0 1
			lset m 2 1 0 1
			lset m 2 3 0 $sq3h
			lset m 3 2 0 $sq3h
		    }
		    Iy {
			lset m 0 1 1 -$sq3h
			lset m 1 0 1 $sq3h
			lset m 1 2 1 -1
			lset m 2 1 1 1
			lset m 2 3 1 -$sq3h
			lset m 3 2 1 $sq3h
		    }
		    Ip {
			lset m 0 1 0 $sq3
			lset m 1 2 0 [expr $s+0.5]
			lset m 2 3 0 $sq3
		    }
		    Im {
			lset m 1 0 0 $sq3
			lset m 2 1 0 [expr $s+0.5]
			lset m 3 2 0 $sq3
			
		    }
		    Ixc {
			lset m 1 2 0 1.0
			lset m 2 1 0 1.0
		    }
		    Iyc {
			lset m 1 2 1 -1.0
			lset m 2 1 1 1.0
		    }
		    Izc {
			lset m 1 1 0 0.5
			lset m 2 2 0 -0.5
		    }
		    Ipc {
			lset m 1 2 0 [expr $s+0.5]
		    }
		    Imc {
			lset m 2 1 0 [expr $s+0.5]
		    }
		    IST1 {
			lset m 0 1 0 $sq3
			lset m 2 3 0 $sq3
		    }
		    IST1p {
			lset m 0 1 0 $sq3
		    }
		    IST1m {
			lset m 2 3 0 $sq3
		    }
		    IST1z {
			lset m 0 0 0 1.5
			lset m 3 3 0 -1.5
		    }
		    IST1x {
			lset m 0 1 0 $sq3h
			lset m 1 0 0 $sq3h
			lset m 2 3 0 $sq3h
			lset m 3 2 0 $sq3h
		    }
		    IST1y {
			lset m 0 1 1 -$sq3h
			lset m 1 0 1 $sq3h
			lset m 2 3 1 -$sq3h
			lset m 3 2 1 $sq3h
		    }
		    I3Qx {
			lset m 0 3 0 0.5
			lset m 3 0 0 0.5
		    }
		    I3Qy {
			lset m 0 3 1 -0.5
			lset m 3 0 1 0.5
		    }
		    I3Qz {
			lset m 0 0 0 0.5
			lset m 3 3 0 -0.5
		    }
		    default {
			puts "Usage: mgetmatrix <spin> <operator>"
			puts "$op ist not a valid operator!"
			puts "Options are: E Ix Iy Iz Ip Im Ixc Iyc Izc Ipc Imc"
			puts "             IST1 IST1p IST1m IST1x IST1y IST1z (satellite transistion)"
			puts "             I3Qx I3Qy I3Qz"
			exit
		    }
		}
	    }
	    5/2 {
		switch -exact -- $op {
		    E {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 1
			}
		    }
		    Iz {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 [expr $s-$i]
			}
		    }
		    Ix {
			lset m 0 1 0 $sq5h
			lset m 1 0 0 $sq5h
			lset m 4 5 0 $sq5h
			lset m 5 4 0 $sq5h
			lset m 1 2 0 $sq2
			lset m 2 1 0 $sq2
			lset m 3 4 0 $sq2
			lset m 4 3 0 $sq2
			lset m 2 3 0 [expr $s-1.0]
			lset m 3 2 0 [expr $s-1.0]
			# lset m row col real/im val
		    }
		    Iy {
			lset m 0 1 1 -$sq5h
			lset m 1 0 1 $sq5h
			lset m 4 5 1 -$sq5h
			lset m 5 4 1 $sq5h
			lset m 1 2 1 -$sq2
			lset m 2 1 1 $sq2
			lset m 3 4 1 -$sq2
			lset m 4 3 1 $sq2
			lset m 2 3 1 [expr -($s-1.0)]
			lset m 3 2 1 [expr $s-1.0]
		    }
		    Ip {
			lset m 0 1 0 $sq5
			lset m 1 2 0 $sq8
			lset m 2 3 0 [expr $s+0.5]
			lset m 3 4 0 $sq8
			lset m 4 5 0 $sq5
		    }
		    Im {
			lset m 1 0 0 $sq5
			lset m 2 1 0 $sq8
			lset m 3 2 0 [expr $s+0.5]
			lset m 4 3 0 $sq8
			lset m 5 4 0 $sq5
		    }
		    Ip {
			lset m 0 1 0 $sq5
			lset m 1 2 0 $sq8
			lset m 2 3 0 [expr $s+0.5]
			lset m 3 4 0 $sq8
			lset m 4 5 0 $sq5
		    }
		    Ipc {
			 lset m 2 3 0 [expr $s+0.5]
		    }
		    Imc {
			lset m 3 2 0 [expr $s+0.5]
		    }
		    Ixc {
			lset m 2 3 0 1.5
			lset m 3 2 0 1.5
		    }
		    Iyc {
			lset m 2 3 1 -1.5
			lset m 3 2 1 1.5
		    }
		    Izc {
			lset m 2 2 0 0.5
			lset m 3 3 0 -0.5
		    }
		    IST1 {
			lset m 1 2 0 $sq8
			lset m 3 4 0 $sq8
		    }
		    IST1p {
			lset m 1 2 0 $sq8
		    }
		    IST1m {
			lset m 3 4 0 $sq8
		    }
		    IST1x {
			lset m 1 2 0 $sq2
			lset m 2 1 0 $sq2
			lset m 3 4 0 $sq2
			lset m 4 3 0 $sq2
		    }
		    IST1y {
			lset m 1 2 1 -$sq2
			lset m 2 1 1 $sq2
			lset m 3 4 1 -$sq2
			lset m 4 3 1 $sq2
		    }
		    IST1z {
			lset m 1 1 0 1.5
			lset m 4 4 0 -1.5
		    }
		    IST2 {
			lset m 0 1 0 $sq5
			lset m 4 5 0 $sq5
		    }
		    IST2p {
			lset m 0 1 0 $sq5
		    }
		    IST2m {
			lset m 4 5 0 $sq5
		    }
		    IST2x {
			lset m 0 1 0 $sq5h
			lset m 1 0 0 $sq5h
			lset m 4 5 0 $sq5h
			lset m 5 4 0 $sq5h
		    }
		    IST2y {
			lset m 0 1 1 -$sq5h
			lset m 1 0 1 $sq5h
			lset m 4 5 1 -$sq5h
			lset m 5 4 1 $sq5h
		    }
		    IST2z {
			lset m 0 0 0 2.5
			lset m 5 5 0 -2.5
		    }
		    default {
			puts "Usage: mgetmatrix <spin> <operator>"
			puts "$op ist not a valid operator!"
			puts "Options are: E Ix Iy Iz Ip Im Ixc Iyc Izc Ipc Imc"
			puts "             IST1 IST1p IST1m IST1x IST1y IST1z (inner satellite transistion)"
			puts "             IST2 IST2p IST2m IST2x IST2y IST2z (outer satellite transistion)"
			#puts "             I3Qx I3Qy I3Qz I5Qx I5Qy I5Qz"
			exit
		    }
		}
	    }
	    7/2 {
		switch -exact -- $op {
		    E {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 1
			}
		    }
		    Iz {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 [expr $s-$i]
			}
		    }
		    Ix {
			lset m 0 1 0 $sq7h
			lset m 1 0 0 $sq7h
			lset m 6 7 0 $sq7h
			lset m 7 6 0 $sq7h
			lset m 1 2 0 $sq3
 			lset m 2 1 0 $sq3
			lset m 5 6 0 $sq3
			lset m 6 5 0 $sq3
			lset m 2 3 0 $sq15h
			lset m 3 2 0 $sq15h
			lset m 4 5 0 $sq15h
			lset m 5 4 0 $sq15h
 			lset m 3 4 0 2
 			lset m 4 3 0 2
		    }
		    Iy {
			lset m 0 1 1 -$sq7h
			lset m 1 0 1 $sq7h
			lset m 6 7 1 -$sq7h
			lset m 7 6 1 $sq7h
			lset m 1 2 1 -$sq3
 			lset m 2 1 1 $sq3
			lset m 5 6 1 -$sq3
			lset m 6 5 1 $sq3
			lset m 2 3 1 -$sq15h
			lset m 3 2 1 $sq15h
			lset m 4 5 1 -$sq15h
			lset m 5 4 1 $sq15h
 			lset m 3 4 1 -2
 			lset m 4 3 1 2

		    }
		    Ip {
			lset m 0 1 0 $sq7
			lset m 1 2 0 $sq12
			lset m 2 3 0 $sq15
			lset m 3 4 0 4
			lset m 4 5 0 $sq15
			lset m 5 6 0 $sq12
			lset m 6 7 0 $sq7
		    }
		    Im {
			lset m 1 0 0 $sq7
			lset m 2 1 0 $sq12
			lset m 3 2 0 $sq15
			lset m 4 3 0 4
			lset m 5 4 0 $sq15
			lset m 6 5 0 $sq12
			lset m 7 6 0 $sq7
		    }
		    Ipc {
			lset m 3 4 0 4.0
		    }
		    Imc {
			lset m 4 3 0 4.0
		    }
		    Ixc {
			lset m 3 4 0 2.0
			lset m 4 3 0 2.0
		    }
		    Iyc {
			lset m 3 4 1 -2.0
			lset m 4 3 1 2.0
		    }
		    Izc {
			lset m 3 3 0 0.5
			lset m 4 4 0 -0.5
		    }
		    IST1 {
			lset m 2 3 0 $sq15
			lset m 4 5 0 $sq15
		    }
		    IST1p {
			lset m 2 3 0 $sq15
		    }
		    IST1m {
			lset m 4 5 0 $sq15
		    }
		    IST1x {
			lset m 2 3 0 $sq15h
			lset m 3 2 0 $sq15h
			lset m 4 5 0 $sq15h
			lset m 5 4 0 $sq15h
		    }
		    IST1y {
			lset m 2 3 1 -$sq15h
			lset m 3 2 1 $sq15h
			lset m 4 5 1 -$sq15h
			lset m 5 4 1 $sq15h
		    }
		    IST1z {
			lset m 2 2 0 1.5
			lset m 5 5 0 -1.5
		    }
		    IST2 {
			lset m 1 2 0 $sq12
			lset m 5 6 0 $sq12
		    }
		    IST2p {
			lset m 1 2 0 $sq12
		    }
		    IST2m {
			lset m 5 6 0 $sq12
		    }
		    IST2x {
			lset m 1 2 0 $sq3
 			lset m 2 1 0 $sq3
			lset m 5 6 0 $sq3
			lset m 6 5 0 $sq3
		    }
		    IST2y {
			lset m 1 2 1 -$sq3
 			lset m 2 1 1 $sq3
			lset m 5 6 1 -$sq3
			lset m 6 5 1 $sq3
		    }
		    IST2z {
			lset m 1 1 0 2.5
			lset m 6 6 0 -2.5
		    }
		    IST3 {
			lset m 1 0 0 $sq7
			lset m 7 6 0 $sq7
		    }
		    IST3p {
			lset m 0 1 0 $sq7
		    }
		    IST3m {
			lset m 6 7 0 $sq7
		    }
		    IST3x {
			lset m 0 1 0 $sq7h
			lset m 1 0 0 $sq7h
			lset m 6 7 0 $sq7h
			lset m 7 6 0 $sq7h
		    }
		    IST3y {
			lset m 0 1 1 -$sq7h
			lset m 1 0 1 $sq7h
			lset m 6 7 1 -$sq7h
			lset m 7 6 1 $sq7h
		    }
		    IST3z {
			lset m 0 0 0 3.5
			lset m 7 7 0 -3.5
		    }
		    default {
			puts "Usage: mgetmatrix <spin> <operator>"
			puts "$op ist not a valid operator!"
			puts "Options are: E Ix Iy Iz Ip Im Ixc Iyc Izc Ipc Imc"
			puts "             IST1 IST1p IST1m IST1x IST1y IST1z (inner satellite transistion)"
			puts "             IST2 IST2p IST2m IST2x IST2y IST2z (middle satellite transistion)"
			puts "             IST3 IST3p IST3m IST3x IST3y IST3z (outer most satellite transistion)"
			exit
		    }
		}   
	    }
	    9/2 {
		switch -exact -- $op {
		    E {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 1
			}
		    }
		    Iz {
			for {set i 0} {$i <= [expr int(2*$s)]} {incr i} {
			    lset m $i $i 0 [expr $s-$i]
			}
		    }
		    Ix {
			lset m 0 1 0 1.5
			lset m 1 0 0 1.5
			lset m 8 9 0 1.5
			lset m 9 8 0 1.5
			lset m 1 2 0 2.0
			lset m 2 1 0 2.0
			lset m 7 8 0 2.0
			lset m 8 7 0 2.0
			lset m 2 3 0 $sq21h
			lset m 3 2 0 $sq21h
			lset m 6 7 0 $sq21h
			lset m 7 6 0 $sq21h
			lset m 3 4 0 $sq6
			lset m 4 3 0 $sq6
			lset m 5 6 0 $sq6
			lset m 6 5 0 $sq6
			lset m 4 5 0 2.5
			lset m 5 4 0 2.5
		    }
		    Iy {
			lset m 0 1 1 -1.5
			lset m 1 0 1 1.5
			lset m 8 9 1 -1.5
			lset m 9 8 1 1.5
			lset m 1 2 1 -2.0
			lset m 2 1 1 2.0
			lset m 7 8 1 -2.0
			lset m 8 7 1 2.0
			lset m 2 3 1 -$sq21h
			lset m 3 2 1 $sq21h
			lset m 6 7 1 -$sq21h
			lset m 7 6 1 $sq21h
			lset m 3 4 1 -$sq6
			lset m 4 3 1 $sq6
			lset m 5 6 1 -$sq6
			lset m 6 5 1 $sq6
			lset m 4 5 1 -2.5
			lset m 5 4 1 2.5
		    }
		    Im {
			lset m 1 0 0 3.0
			lset m 2 1 0 4.0
			lset m 3 2 0 $sq21
			lset m 4 3 0 $sq24
			lset m 5 4 0 5.0
			lset m 6 5 0 $sq24
			lset m 7 6 0 $sq21
			lset m 8 7 0 4.0
			lset m 9 8 0 3.0
		    }
		    Ip {
			lset m 0 1 0 3.0
			lset m 1 2 0 4.0
			lset m 2 3 0 $sq21
			lset m 3 4 0 $sq24
			lset m 4 5 0 5.0
			lset m 5 6 0 $sq24
			lset m 6 7 0 $sq21
			lset m 7 8 0 4.0
			lset m 8 9 0 3.0
		    }
		    Ipc {
			lset m 4 5 0 5.0
		    }
		    Imc {
			lset m 5 4 0 5.0
		    }
		    Ixc {
			lset m 4 5 0 2.5
			lset m 5 4 0 2.5
		    }
		    Iyc {
			lset m 4 5 1 -2.5
			lset m 5 4 1 2.5
		    }
		    Izc {
			lset m 4 4 0 0.5
			lset m 5 5 0 -0.5
		    }
		    IST1 {
			lset m 3 4 0 $sq24
			lset m 5 6 0 $sq24
		    }
		    IST1p {
			lset m 3 4 0 $sq24			
		    }
		    IST1m {
			lset m 5 6 0 $sq24
		    }
		    IST1x {
			lset m 3 4 0 $sq6
			lset m 4 3 0 $sq6
			lset m 5 6 0 $sq6
			lset m 6 5 0 $sq6
		    }
		    IST1y {
			lset m 3 4 1 -$sq6
			lset m 4 3 1 $sq6
			lset m 5 6 1 -$sq6
			lset m 6 5 1 $sq6
		    }
		    IST1z {
			lset m 3 3 0 1.5
			lset m 6 6 0 -1.5
		    }
		    IST2 {
			lset m 2 3 0 $sq21
			lset m 6 7 0 $sq21
		    }
		    IST2p {
			lset m 2 3 0 $sq21
		    }
		    IST2m {
			lset m 6 7 0 $sq21
		    }
		    IST2x {
			lset m 2 3 0 $sq21h
			lset m 3 2 0 $sq21h
			lset m 6 7 0 $sq21h
			lset m 7 6 0 $sq21h
		    }
		    IST2y {
			lset m 2 3 1 -$sq21h
			lset m 3 2 1 $sq21h
			lset m 6 7 1 -$sq21h
			lset m 7 6 1 $sq21h
		    }
		    IST2z {
			lset m 2 2 0 2.5
			lset m 7 7 0 -2.5
		    }
		    IST3 {
			lset m 1 2 0 4.0
			lset m 7 8 0 4.0
		    }
		    IST3p {
			lset m 1 2 0 4.0		
		    }
		    IST3m {
			lset m 7 8 0 4.0
		    }
		    IST3x {
			lset m 1 2 0 2.0
			lset m 2 1 0 2.0
			lset m 7 8 0 2.0
			lset m 8 7 0 2.0
		    }
		    IST3y {
			lset m 1 2 1 -2.0
			lset m 2 1 1 2.0
			lset m 7 8 1 -2.0
			lset m 8 7 1 2.0
		    }
		    IST3z {
			lset m 1 1 0 3.5
			lset m 8 8 0 -3.5
		    }
		    IST4 {
			lset m 0 1 0 3.0
			lset m 8 9 0 3.0
		    }
		    IST4p {
			lset m 0 1 0 3.0
		    }
		    IST4m {
			lset m 8 9 0 3.0
		    }
		    IST4x {
			lset m 0 1 0 1.5
			lset m 1 0 0 1.5
			lset m 8 9 0 1.5
			lset m 9 8 0 1.5
		    }
		    IST4y {
			lset m 0 1 1 -1.5
			lset m 1 0 1 1.5
			lset m 8 9 1 -1.5
			lset m 9 8 1 1.5
		    }
		    IST4z {
			lset m 0 0 0 4.5
			lset m 9 9 0 -4.5
		    }
		    default {
			puts "Usage: mgetmatrix <spin> <operator>"
			puts "$op ist not a valid operator!"
			puts "Options are: E Ix Iy Iz Ip Im Ixc Iyc Izc Ipc Imc"
			puts "             IST1 IST1p IST1m IST1x IST1y IST1z (inner most satellite transistion)"
			puts "             IST2 IST2p IST2m IST2x IST2y IST2z (satellite transistion)"
			puts "             IST3 IST3p IST3m IST3x IST3y IST3z (satellite transistion)"
			puts "             IST4 IST4p IST4m IST4x IST4y IST4z (outer most satellite transistion)"
			exit
		    }
		}
	    }
	}
	return $m
    }

    proc mdirprod {m1 m2} {
	set r1 [llength $m1]
	set c1 [llength [lindex $m1 0]]
	set r2 [llength $m2]
	set c2 [llength [lindex $m2 0]]
	set rows [expr $r1*$r2]
	set cols [expr $c1*$c2]
	set m [mcreatematrix $rows $cols]
	for {set r 0} {$r < $rows} {incr r} {
	    for {set c 0} {$c < $cols} {incr c} {
		set m1re [lindex $m1 [expr ($r/$r2)] [expr int($c/$c2)] 0]
		set m2re [lindex $m2 [expr $r % $r1] [expr $c % $c1] 0]
		set m1im [lindex $m1 [expr ($r/$r2)] [expr int($c/$c2)] 1]
		set m2im [lindex $m2 [expr $r % $r1] [expr $c % $c1] 1]
                lset m $r $c 0 [expr $m1re*$m2re-$m1im*$m2im]
                lset m $r $c 1 [expr $m1re*$m2im+$m1im*$m2re]
		#lset m $r $c 0 [expr $m1re*$m2re]
		#lset m $r $c 1 [expr $m1im*$m2im]
	    }
	}
	return $m
    }
}

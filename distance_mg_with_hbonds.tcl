# calculating the out of plane distance of mg in the final structure

# matrix order

proc mg {} {

	# GETTING THE COORDINATES FROM THE TRAJECTORY FILE TO GET THE GOOD ENSEMBLE OF VALUES

	set f [open "frame.mdcrd" "r"]
	set data [read $f]
	close $f

	set num_cla [input2]
	set num_per_res [input3]

	# MG COORDINATES

	set atom_pos 0
	set elem ""

	set start [expr { 3 + $atom_pos }]
	set end [expr { $start + (3*$num_cla*$num_per_res) }]
	set count $start
	set res 1
	for {set i $start} {$i < $end} {incr i [expr { 3*$num_per_res }]} {
		set coordx [lindex $data $i]
		set coordy [lindex $data [expr { $i + 1 }]]
		set coordz [lindex $data [expr { $i + 2 }]]
		set elem2 [list $coordx $coordy $coordz $res]
		set elem [linsert $elem end $elem2]
		incr res
	}
	return $elem
}

proc com {} {

		# GETTING THE COORDINATES FROM THE TRAJECTORY FILE TO GET THE GOOD ENSEMBLE OF VALUES

	set f [open "frame.mdcrd" "r"]
	set data [read $f]
	close $f

	set g [open "[input4]" "r"]
	set data1 [read $g]
	close $g

	set h [open "frame.pdb" "r"]
	set data2 [read $h]
	close $h

	set k1 0

	# READING THE MASS FROM THE PRMTOP FILE

	while { [lindex $data1 $k1] != "MASS" } {
		incr k1
	}
	incr k1 2
	set k2 $k1

	set num_cla [input2]
	set num_per_res [input3]
	set num_frames 100

	set nres_atoms [expr { $num_cla * $num_per_res }]

	# CENTRE OF MASS COORDINATES

	set atom_pos 0
	set elem ""

	set start [expr { 3 + $atom_pos }]
	set end [expr { $start + (3*$num_cla*$num_per_res) }]
	set k1 $k2
	set comx 0.0
	set comy 0.0
	set comz 0.0
	set total_mass 0.0
	set count $start
	set res 1
	for {set i $start} {$i < $end} {incr i 3} {
		set mass [lindex $data1 $k1]
		set total_mass [expr {$total_mass + $mass}]
		incr k1
		set coordx [lindex $data $i]
		set coordy [lindex $data [expr { $i + 1 }]]
		set coordz [lindex $data [expr { $i + 2 }]]
		set comx [expr { $comx + ($mass*$coordx) }]
		set comy [expr { $comy + ($mass*$coordy) }]
		set comz [expr { $comz + ($mass*$coordz) }]
		if { [expr { $count % (3*$num_per_res) }] == 0 } {
			set comx [expr { $comx / $total_mass }]
			set comy [expr { $comy / $total_mass }]
			set comz [expr { $comz / $total_mass }]
			set elem2 [list $comx $comy $comz $res]
			set elem [linsert $elem end $elem2]
			set comx 0.0
			set comy 0.0
			set comz 0.0
			set total_mass 0.0
			set count $start
			incr res
			set k1 $k2
		}
		incr count 3		
	}
	return $elem
}
	
proc mdcrd_coord {atom_pos} {
	
	# GETTING THE COORDINATES FROM THE TRAJECTORY FILE TO GET THE GOOD ENSEMBLE OF VALUES

	set f [open "frame.mdcrd" "r"]
	set data [read $f]
	close $f

	set g [open "[input4]" "r"]
	set data1 [read $g]
	close $g

	set h [open "frame.pdb" "r"]
	set data2 [read $h]
	close $h

	set h1 [open "dummy" "w"]

	# GETTING THE POSITION OF THE NITROGENS WITH RESPECT TO MG

	set k 0

	while {[lindex $data2 $k] != "TER" || [lindex $data2 [expr { $k - 1 }]] != "H" } {
		if { [lindex $data2 $k] == "N" } {
			set pos [lindex $data2 [expr { $k - 9 }]]
			puts $h1 "$pos"
		}
	incr k
	}
	close $h1 
	set h1 [open "dummy" "r"]
	set dum [read $h1]
	close $h1

	set n1 [expr { [lindex $dum 0] -1 }] 
	set n2 [expr { [lindex $dum 1] -1 }]
	set n3 [expr { [lindex $dum 2] -1 }]
	set n4 [expr { [lindex $dum 3] -1 }]

	set k1 0

	# READING THE MASS FROM THE PRMTOP FILE

	while { [lindex $data1 $k1] != "MASS" } {
		incr k1
	}
	incr k1 2

	set num_cla [input2]
	set num_per_res [input3]
	set num_frames 100

	set nres_atoms [expr { $num_cla * $num_per_res }]

	# CENTRE OF MASS COORDINATES

	set comx 0.0
	set comy 0.0
	set comz 0.0
	set total_mass 0.0

	set start [expr { 3 + $atom_pos }]
	set end [expr { $start + (3*$num_per_res) }]
	for {set i $start} {$i < $end} {incr i 3} {
			set mass [lindex $data1 $k1]
			set total_mass [expr {$total_mass + $mass}]
			incr k1
			set coordx [lindex $data $i]
			set coordy [lindex $data [expr { $i + 1 }]]
			set coordz [lindex $data [expr { $i + 2 }]]
			set comx [expr { $comx + ($mass*$coordx) }]
			set comy [expr { $comy + ($mass*$coordy) }]
			set comz [expr { $comz + ($mass*$coordz) }]
		}
	set comx [expr { $comx / $total_mass }]
	set comy [expr { $comy / $total_mass }]
	set comz [expr { $comz / $total_mass }] 

	set n 0
	for {set i $start} {$i < $end} {incr i 3} {
		if { $i == $start } {
			set cord_mgx [expr { [lindex $data $i] - $comx }]
			set cord_mgy [expr { [lindex $data [expr { $i + 1 }]] - $comy }]
			set cord_mgz [expr { [lindex $data [expr { $i + 2 }]] - $comz }]
			
		}
		if { $i == [expr { $start + (3*$n1) }] || $i == [expr { $start + (3*$n2) }] || $i == [expr { $start + (3*$n3) }] || $i == [expr { $start + (3*$n4) }]} {
			set cord_x($n) [expr { [lindex $data $i] - $comx }]
			set cord_y($n) [expr { [lindex $data [expr { $i + 1}]] - $comy }]
			set cord_z($n) [expr { [lindex $data [expr { $i + 2 }]] - $comz }]
			incr n
		}

		# Phytol ester

		if { $i == [expr { $start + 42 }] } {
			set normx [expr { [lindex $data [expr { $i + 3 }]] - $comx }] 
			set normy [expr { [lindex $data [expr { $i + 4 }]] - $comy }] 
			set normz [expr { [lindex $data [expr { $i + 5 }]] - $comz }] 
		}
	}
	#set normx [expr { [lindex $data [expr { $end - 3 }]] - $comx }]
	#set normy [expr { [lindex $data [expr { $end - 2 }]] - $comy }]
	#set normz [expr { [lindex $data [expr { $end - 1 }]] - $comz }]

	return [list $cord_mgx $cord_mgy $cord_mgz $cord_x(0) $cord_y(0) $cord_z(0) $cord_x(1) $cord_y(1) $cord_z(1) $cord_x(2) $cord_y(2) $cord_z(2) $cord_x(3) $cord_y(3) $cord_z(3) $normx $normy $normz]
}	

proc elem_cord {a} {

	# GETTING THE COORDINATES WITH RESPECT TO CENTRE OF MASS EACH RESIDUE

	set f [open "frame.pdb" "r"]
	set data [read $f]
	close $f

	set g [open "[input4]" "r"]
	set data1 [read $g]
	close $g

	set k $a
	set k2 $a

	set i 0
	set t 0

	set k1 0

	# READING THE MASS FROM THE PRMTOP FILE

	while { [lindex $data1 $k1] != "MASS" } {
		incr k1
	}
	incr k1 2

	# CENTRE OF MASS COORDINATES

	set comx 0.0
	set comy 0.0
	set comz 0.0
	set total_mass 0.0

	while {[lindex $data $k2] != "TER" || [lindex $data [expr { $k2 - 1 }]] != "H" } {
		if { [lindex $data $k2] == "ATOM" } { 
			set mass [lindex $data1 $k1]
			set total_mass [expr {$total_mass + $mass}]
			incr k1
			set coordx [lindex $data [expr { $k2 + 5 }]] 
			set coordy [lindex $data [expr { $k2 + 6 }]] 
			set coordz [lindex $data [expr { $k2 + 7 }]]
			set comx [expr { $comx + ($mass*$coordx) }]
			set comy [expr { $comy + ($mass*$coordy) }]
			set comz [expr { $comz + ($mass*$coordz) }]
		}
		incr k2
	}
	set comx [expr { $comx / $total_mass }]
	set comy [expr { $comy / $total_mass }]
	set comz [expr { $comz / $total_mass }] 
	
	while {[lindex $data $k] != "TER" || [lindex $data [expr { $k - 1 }]] != "H" } {
		if { [lindex $data $k] == "MG" } {
			set cord_mgx [expr { [lindex $data [expr { $k - 5 }]] - $comx }]
			set cord_mgy [expr { [lindex $data [expr { $k - 4 }]] - $comy }]
			set cord_mgz [expr { [lindex $data [expr { $k - 3 }]] - $comz }]
		}
		if { [lindex $data $k] == "N" } {
			set cord_x($i) [expr { [lindex $data [expr { $k - 5 }]] - $comx }]
			set cord_y($i) [expr { [lindex $data [expr { $k - 4 }]] - $comy }]
			set cord_z($i) [expr { [lindex $data [expr { $k - 3 }]] - $comz }]
			incr i
		}
		incr k
	}
	set normx [expr { [lindex $data [expr { $k - 6 }]] - $comx }]
	set normy [expr { [lindex $data [expr { $k - 5 }]] - $comy }]
	set normz [expr { [lindex $data [expr { $k - 4 }]] - $comz }]

	return [list $cord_mgx $cord_mgy $cord_mgz $cord_x(0) $cord_y(0) $cord_z(0) $cord_x(1) $cord_y(1) $cord_z(1) $cord_x(2) $cord_y(2) $cord_z(2) $cord_x(3) $cord_y(3) $cord_z(3) $k $normx $normy $normz]
}

proc plane_N {} {

	package require math::linearalgebra 

	set n_cla [input2]

	set num [open "[input4]" "r"]	
	set natom [read $num]
	close $num

	set k 0
	while { [lindex $natom $k] != "%FORMAT(10I8)" } {
		incr k
	}
	incr k
	set n_atoms [lindex $natom $k]
	set num_at_pres [input3]

	for {set i 0} {$i < $n_cla} {incr i} {
		set dis($i) 0.0
	}
	
	set m [open "[output2]" "w"]
	set pa [open "angle_plane" "w"]

	set start_frame [input5]
	set end_frame [input6]
	set step [input7]
	set n_frames [expr { (($end_frame - $start_frame)/$step) + 1 }]

	for {set fr $start_frame} {$fr <= $end_frame} {incr fr $step} {

		set pos 0

		# EXECUTING CPPTRAJ

		set inp [open "ai.dat" "w"]
		puts $inp "trajin [cpptraj_traj] $fr $fr"
		puts $inp "trajout frame.mdcrd mdcrd nobox"
		puts $inp "go"
		close $inp

		set inp [open "ai.dat" "r"]

		exec cpptraj -p [input4] -i ai.dat 

		close $inp

		for {set i 0} {$i < $n_cla} {incr i} {
			puts "		****** FRAME $fr :: step [expr { $i + 1 }] of $n_cla :: of [input6] frames   ******		"
			if { $i == 0 } {
				set f [open "$i.avg" "w"]
				set coord [mdcrd_coord $pos]
				puts $f "{$coord}"
				close $f
			} else {
				set pos [expr { $pos + (3*$num_at_pres) }]
				set f [open "$i.avg" "w"]
				set coord [mdcrd_coord $pos]
				puts $f "{$coord}"
				close $f
			}
		}

		# determining the equation of plane taking three nitrogens at a time ax+by+cz+d=0 for all CLA and the distance of this plane from mg the mean of this plane

		# A matrix
		for {set i 0} {$i < $n_cla} {incr i} {
			#puts "		***** Determing distance of mg in residue [expr {$i + 1}] ***** 		"
			set f [open "$i.avg" "r"]
			set data [read $f]
			close $f

			# PHYTOL HYDROGENS

			set phy_H [list [lindex $data 0 15] [lindex $data 0 16] [lindex $data 0 17]]

			# vector 1-2

			set elem_00 [expr { [lindex $data 0 6] - [lindex $data 0 3] }]
			set elem_01 [expr { [lindex $data 0 7] - [lindex $data 0 4] }]
			set elem_02 [expr { [lindex $data 0 8] - [lindex $data 0 5] }]

			# vector 2-3

			set elem_10 [expr { [lindex $data 0 9] - [lindex $data 0 6] }]
			set elem_11 [expr { [lindex $data 0 10] - [lindex $data 0 7] }]
			set elem_12 [expr { [lindex $data 0 11] - [lindex $data 0 8] }]

			# vector 3-1

			set elem_20 [expr { [lindex $data 0 3] - [lindex $data 0 9] }]
			set elem_21 [expr { [lindex $data 0 4] - [lindex $data 0 10] }]
			set elem_22 [expr { [lindex $data 0 5] - [lindex $data 0 11] }]

			# vector 4-1

			set elem_30 [expr { [lindex $data 0 3] - [lindex $data 0 12] }]
			set elem_31 [expr { [lindex $data 0 4] - [lindex $data 0 13] }]
			set elem_32 [expr { [lindex $data 0 5] - [lindex $data 0 14] }]

			# vector 4-2

			set elem_40 [expr { [lindex $data 0 6] - [lindex $data 0 12] }]
			set elem_41 [expr { [lindex $data 0 7] - [lindex $data 0 13] }]
			set elem_42 [expr { [lindex $data 0 8] - [lindex $data 0 14] }]

			# vector 3-4

			set elem_50 [expr { [lindex $data 0 12] - [lindex $data 0 9] }]
			set elem_51 [expr { [lindex $data 0 13] - [lindex $data 0 10] }]
			set elem_52 [expr { [lindex $data 0 14] - [lindex $data 0 11] }]
		
			# solution 1 2 3

			set vec1 [list $elem_00 $elem_01 $elem_02]
			set vec2 [list $elem_10 $elem_11 $elem_12]
			set vec3 [list [lindex $data 0 3] [lindex $data 0 4] [lindex $data 0 5]]
			set norm_vect1 [::math::linearalgebra::crossproduct $vec2 $vec1]

			if { [::math::linearalgebra::dotproduct $norm_vect1 $phy_H] < 0.0 } {
				set norm_vect1 [::math::linearalgebra::scale -1.0 $norm_vect1]
			}
			set norm_vect1 [::math::linearalgebra::unitLengthVector $norm_vect1]
			set d1 [::math::linearalgebra::dotproduct $norm_vect1 $vec3]

			# solution 2 3 4

			set vec1 [list $elem_10 $elem_11 $elem_12]
			set vec2 [list $elem_50 $elem_51 $elem_52]
			set vec3 [list [lindex $data 0 6] [lindex $data 0 7] [lindex $data 0 8]]
			set norm_vect2 [::math::linearalgebra::crossproduct $vec2 $vec1]

			if { [::math::linearalgebra::dotproduct $norm_vect2 $phy_H] < 0.0 } {
				set norm_vect2 [::math::linearalgebra::scale -1.0 $norm_vect2]
			}
			set norm_vect2 [::math::linearalgebra::unitLengthVector $norm_vect2]
			set d2 [::math::linearalgebra::dotproduct $norm_vect2 $vec3]

			# solution 3 4 1

			set vec1 [list $elem_50 $elem_51 $elem_52]
			set vec2 [list $elem_30 $elem_31 $elem_32]
			set vec3 [list [lindex $data 0 3] [lindex $data 0 4] [lindex $data 0 5]]
			set norm_vect3 [::math::linearalgebra::crossproduct $vec2 $vec1]

			if { [::math::linearalgebra::dotproduct $norm_vect3 $phy_H] < 0.0 } {
				set norm_vect3 [::math::linearalgebra::scale -1.0 $norm_vect3]
			}
			set norm_vect3 [::math::linearalgebra::unitLengthVector $norm_vect3]
			set d3 [::math::linearalgebra::dotproduct $norm_vect3 $vec3]

			# solution 4 1 2

			set vec1 [list $elem_30 $elem_31 $elem_32]
			set vec2 [list $elem_00 $elem_01 $elem_02]
			set vec3 [list [lindex $data 0 3] [lindex $data 0 4] [lindex $data 0 5]]
			set norm_vect4 [::math::linearalgebra::crossproduct $vec2 $vec1]

			if { [::math::linearalgebra::dotproduct $norm_vect4 $phy_H] < 0.0 } {
				set norm_vect4 [::math::linearalgebra::scale -1.0 $norm_vect4]
			}
			set norm_vect4 [::math::linearalgebra::unitLengthVector $norm_vect4]
			set d4 [::math::linearalgebra::dotproduct $norm_vect4 $vec3]

			# determing the equation of the mean plane contaiing all 4 nitrogens

			set norm_vect [::math::linearalgebra::add $norm_vect1 $norm_vect2]
			set norm_vect [::math::linearalgebra::add $norm_vect $norm_vect3]
			set norm_vect [::math::linearalgebra::add $norm_vect $norm_vect4]

			set norm_vect [::math::linearalgebra::scale 0.25 $norm_vect]

			set norm_vect [::math::linearalgebra::unitLengthVector $norm_vect]

			set d [expr { ($d1 + $d2 + $d3 + $d4) / 4.0 }]

			if { $i == 22 || $i == 23 } {
				set elem($i) $norm_vect
			}
	 
			set h [open "dummy" "w"]
			puts $h "$norm_vect"
			puts $h "$d"
			close $h

			set h [open "dummy" "r"]
			set data3 [read $h]
			close $h

			set a [lindex $data3 0]
			set a2 [expr { $a * $a }]
			set b [lindex $data3 1]
			set b2 [expr { $b * $b }]
			set c [lindex $data3 2]
			set c2 [expr { $c * $c }]
			set d [expr { [lindex $data3 3] * -1 }]

			# determing the distance of mg from the chlorophyll center

			set mg [list [lindex $data 0 0] [lindex $data 0 1] [lindex $data 0 2]]

			set mg_pos [::math::linearalgebra::dotproduct $mg $phy_H]
			if { $mg_pos < 0.0 } {
				set position 0
			} else {
				set position 1
			}
	
			set num [::math::linearalgebra::dotproduct $norm_vect $mg]
	
			set num [expr {$num + $d}]
			set num [expr { abs($num) }]
			set den [expr { sqrt($a2 + $b2 + $c2) }]

			set dis($i) [expr { $dis($i) + ($num / $den) }]
		}
		set ang [::math::linearalgebra::dotproduct $elem(22) $elem(23)]
		set ang [expr { acos($ang) }]
		set ang [expr { ($ang * 180) / 3.14 }]
		puts $pa "$fr $ang"
	}

	for {set i 0} {$i < $n_cla} {incr i} {
		set avg_dis [expr { $dis($i) / $n_frames }]
		puts $m "$i $avg_dis"
	}
	close $m
	close $pa

}

proc plane_N_static {} {

	package require math::linearalgebra 

	# EXECUTING CPPTRAJ

	set inp [open "ai.dat" "w"]
	puts $inp "trajin [cpptraj_traj] [input6] [input6]"
	puts $inp "trajout frame.pdb pdb"
	puts $inp "go"
	close $inp

	set inp [open "ai.dat" "r"]

	exec cpptraj -p [input4] -i ai.dat 

	close $inp

	set n_cla [input2]

	for {set i 0} {$i < $n_cla} {incr i} {
		puts "		****** step [expr { $i + 1 }] of $n_cla ******		"
		if { $i == 0 } {
			set f [open "$i" "w"]
			set coord [elem_cord 0]
			puts $f "{$coord}"
			close $f
		} else {
			set g [open "[expr { $i - 1 }]" "r"]
			set data [read $g]
			close $g
			set f [open "$i" "w"]
			set k [expr { [lindex $data 0 15] + 4 }]
			set coord [elem_cord $k]
			puts $f "{$coord}"
			close $f
		}
	}

	# determining the equation of plane taking three nitrogens at a time ax+by+cz+d=0 for all CLA and the distance of this plane from mg the mean of this plane

	# A matrix
	set m [open "[output1]" "w"]
	for {set i 0} {$i < $n_cla} {incr i} {
		puts "		***** Determing distance of mg in residue [expr {$i + 1}] ***** 		"
		set f [open "$i" "r"]
		set data [read $f]
		close $f

		# PHYTOL HYDROGENS

		set phy_H [list [lindex $data 0 16] [lindex $data 0 17] [lindex $data 0 18]]

		# vector 1-2

		set elem_00 [expr { [lindex $data 0 6] - [lindex $data 0 3] }]
		set elem_01 [expr { [lindex $data 0 7] - [lindex $data 0 4] }]
		set elem_02 [expr { [lindex $data 0 8] - [lindex $data 0 5] }]

		# vector 2-3

		set elem_10 [expr { [lindex $data 0 9] - [lindex $data 0 6] }]
		set elem_11 [expr { [lindex $data 0 10] - [lindex $data 0 7] }]
		set elem_12 [expr { [lindex $data 0 11] - [lindex $data 0 8] }]

		# vector 3-1

		set elem_20 [expr { [lindex $data 0 3] - [lindex $data 0 9] }]
		set elem_21 [expr { [lindex $data 0 4] - [lindex $data 0 10] }]
		set elem_22 [expr { [lindex $data 0 5] - [lindex $data 0 11] }]

		# vector 4-1

		set elem_30 [expr { [lindex $data 0 3] - [lindex $data 0 12] }]
		set elem_31 [expr { [lindex $data 0 4] - [lindex $data 0 13] }]
		set elem_32 [expr { [lindex $data 0 5] - [lindex $data 0 14] }]

		# vector 4-2

		set elem_40 [expr { [lindex $data 0 6] - [lindex $data 0 12] }]
		set elem_41 [expr { [lindex $data 0 7] - [lindex $data 0 13] }]
		set elem_42 [expr { [lindex $data 0 8] - [lindex $data 0 14] }]

		# vector 3-4

		set elem_50 [expr { [lindex $data 0 12] - [lindex $data 0 9] }]
		set elem_51 [expr { [lindex $data 0 13] - [lindex $data 0 10] }]
		set elem_52 [expr { [lindex $data 0 14] - [lindex $data 0 11] }]
		
		# solution 1 2 3

		set vec1 [list $elem_00 $elem_01 $elem_02]
		set vec2 [list $elem_10 $elem_11 $elem_12]
		set vec3 [list [lindex $data 0 3] [lindex $data 0 4] [lindex $data 0 5]]
		set norm_vect1 [::math::linearalgebra::crossproduct $vec2 $vec1]

		if { [::math::linearalgebra::dotproduct $norm_vect1 $phy_H] < 0.0 } {
			set norm_vect1 [::math::linearalgebra::scale -1.0 $norm_vect1]
		}

		set norm_vect1 [::math::linearalgebra::unitLengthVector $norm_vect1]
				
		set d1 [::math::linearalgebra::dotproduct $norm_vect1 $vec3]

		# solution 2 3 4

		set vec1 [list $elem_10 $elem_11 $elem_12]
		set vec2 [list $elem_50 $elem_51 $elem_52]
		set vec3 [list [lindex $data 0 6] [lindex $data 0 7] [lindex $data 0 8]]
		set norm_vect2 [::math::linearalgebra::crossproduct $vec2 $vec1]

		if { [::math::linearalgebra::dotproduct $norm_vect2 $phy_H] < 0.0 } {
			set norm_vect2 [::math::linearalgebra::scale -1.0 $norm_vect2]
		}
		set norm_vect2 [::math::linearalgebra::unitLengthVector $norm_vect2]
		set d2 [::math::linearalgebra::dotproduct $norm_vect2 $vec3]

		# solution 3 4 1

		set vec1 [list $elem_50 $elem_51 $elem_52]
		set vec2 [list $elem_30 $elem_31 $elem_32]
		set vec3 [list [lindex $data 0 3] [lindex $data 0 4] [lindex $data 0 5]]
		set norm_vect3 [::math::linearalgebra::crossproduct $vec2 $vec1]

		if { [::math::linearalgebra::dotproduct $norm_vect3 $phy_H] < 0.0 } {
			set norm_vect3 [::math::linearalgebra::scale -1.0 $norm_vect3]
		}
		set norm_vect3 [::math::linearalgebra::unitLengthVector $norm_vect3]
		
		set d3 [::math::linearalgebra::dotproduct $norm_vect3 $vec3]

		# solution 4 1 2

		set vec1 [list $elem_30 $elem_31 $elem_32]
		set vec2 [list $elem_00 $elem_01 $elem_02]
		set vec3 [list [lindex $data 0 3] [lindex $data 0 4] [lindex $data 0 5]]
		set norm_vect4 [::math::linearalgebra::crossproduct $vec2 $vec1]

		if { [::math::linearalgebra::dotproduct $norm_vect4 $phy_H] < 0.0 } {
			set norm_vect4 [::math::linearalgebra::scale -1.0 $norm_vect4]
		}
		set norm_vect4 [::math::linearalgebra::unitLengthVector $norm_vect4]
		set d4 [::math::linearalgebra::dotproduct $norm_vect4 $vec3]

		# determing the equation of the mean plane contaiing all 4 nitrogens

		set norm_vect [::math::linearalgebra::add $norm_vect1 $norm_vect2]
		set norm_vect [::math::linearalgebra::add $norm_vect $norm_vect3]
		set norm_vect [::math::linearalgebra::add $norm_vect $norm_vect4]

		set norm_vect [::math::linearalgebra::scale 0.25 $norm_vect]

		set norm_vect [::math::linearalgebra::unitLengthVector $norm_vect]

		if { $i == 3 || $i == 10 } {
			set elem($i) $norm_vect
		}

		set d [expr { ($d1 + $d2 + $d3 + $d4) / 4.0 }]
 
		set h [open "dummy" "w"]
		puts $h "$norm_vect"
		puts $h "$d"
		close $h

		set h [open "dummy" "r"]
		set data3 [read $h]
		close $h

		set a [lindex $data3 0]
		set a2 [expr { $a * $a }]
		set b [lindex $data3 1]
		set b2 [expr { $b * $b }]
		set c [lindex $data3 2]
		set c2 [expr { $c * $c }]
		set d [expr { [lindex $data3 3] * -1 }]

		# determing the distance of mg from the chlorophyll center

		set mg [list [lindex $data 0 0] [lindex $data 0 1] [lindex $data 0 2]]

		set mg_pos [::math::linearalgebra::dotproduct $mg $phy_H]
		if { $mg_pos < 0.0 } {
			set position 0
		} else {
			set position 1
		}
	
		set num [::math::linearalgebra::dotproduct $norm_vect $mg]
	
		set num [expr {$num + $d}]
		set num [expr { abs($num) }]
		set den [expr { sqrt($a2 + $b2 + $c2) }]

		set dis [expr { $num / $den }]

		puts $m "$i $dis $position"
	}
	close $m
	puts "[::math::linearalgebra::dotproduct $elem(3) $elem(10)]"
}

# THE PROCEDURES AFTER THIS WILL LOOK FOR THE WATER WITHIN 2.5A FROM WACH MG AND CALCULATE THE AVERAGE DISTANCE OF MG FROM WATER, AVERAGED OVER THE EQUILIBRIUM CONFIGURATION

proc water_distance {} {
	
	# DETERMINING THE NEAREST WATER NEIGBOUR TO THE MG IN EACH RESIDUES

	set n_cla [input2]	

	for {set i 1} {$i <= $n_cla} {incr i} {
		set d($i) 0.0
		set ep1($i) 0.0
		set ep2($i) 0.0
		set wat_l($i) ""
	}
	
	set g [open "water_dis" "w"]
	set g1 [open "ester_dis" "w"]

	set start_frame [input5]
	set end_frame [input6]
	set step [input7]
	set n_frames [expr { (($end_frame - $start_frame)/$step) + 1 }]
	for {set fr $start_frame} {$fr <= $end_frame} {incr fr $step} {

		# EXECUTING CPPTRAJ

		set inp [open "ai.dat" "w"]
		puts $inp "trajin [cpptraj_traj] $fr $fr"
		puts $inp "trajout frame.pdb pdb nobox"
		puts $inp "trajout frame.rst rst"
		puts $inp "go"
		close $inp

		set inp [open "ai.dat" "r"]

		exec cpptraj -p [input4] -i ai.dat 

		close $inp
		
		set f1 [open "frame.rst" "r"]
		set data1 [read $f1]
		close $f1

		set box_x [lindex $data1 [expr { [llength $data1] - 3 }]]
		set box_y [lindex $data1 [expr { [llength $data1] - 2 }]]
		set box_z [lindex $data1 [expr { [llength $data1] - 1 }]]

		set f [open "frame.pdb" "r"]
		set data [read $f]
		close $f

		for {set j 1} {$j <= $n_cla} {incr j} {

			set d_temp 0.0

			puts "		****** FRAME $fr :: step $j of $n_cla :: of [input6] frames   ******		"

			set k 0
			while { [lindex $data $k] != $j || [lindex $data [expr { $k - 2 }]] != "Mg1" } {
				incr k
			}
			set mgx [lindex $data [expr { $k + 1 }]]
			set mgy [lindex $data [expr { $k + 2 }]]
			set mgz [lindex $data [expr { $k + 3 }]]

			set k 0

			set temp_dis1 1000.0
			set temp_dis2 1000.0
			set pair1 1000000

			while { [lindex $data $k] != "END" } {
				if { [lindex $data $k] == "O" && [lindex $data [expr { $k + 1 }]] == "WAT" } {
					set ox [lindex $data [expr { $k + 3 }]]
					set oy [lindex $data [expr { $k + 4 }]]
					set oz [lindex $data [expr { $k + 5 }]]

					set difx [expr { $mgx-$ox }]
					if { $difx > [expr { $box_x / 2.0 }] } {
						set difx [expr { $difx - $box_x }]
					}
					set difx [expr { $difx * $difx }]

					set dify [expr { $mgy-$oy }]
					if { $dify > [expr { $box_y / 2.0 }] } {
						set dify [expr { $dify - $box_y }]
					}
					set dify [expr { $dify * $dify }]

					set difz [expr { $mgz-$oz }]
					if { $difz > [expr { $box_z / 2.0 }] } {
						set difz [expr { $difz - $box_z }]
					}
					set difz [expr { $difz * $difz }]

					set dis [expr { sqrt($difx + $dify + $difz) }]
		
					if { $dis < $temp_dis1 } {
							set temp_dis2 $temp_dis1
							set temp_dis1 $dis
							set pair2 $pair1
							set pair1 $k					
					} elseif { $dis < $temp_dis2 } {
							set temp_dis2 $dis
							set pair2 $k
					} 
				}
				incr k 
			}

			# CALCULATING THE DISTANCE OF THE METHYL ESTER FROM THE HYDROGENS OF THE BONDED WATER MOLECULE

			ester_dis $j [expr { $pair1 - 1 }] [expr { $pair2 - 1 }]

			set dum_es [open "dummy_ester" "r"]
			set de [read $dum_es]
			close $dum_es

			set ep1($j) [expr { $ep1($j) + [lindex $de 0] }]
			set ep2($j) [expr { $ep2($j) + [lindex $de 1] }]

			if { [expr { $fr % 10 }] == 0 } {
				set wat_l($j) [linsert $wat_l($j) end [lindex $data [expr { $pair1 + 2 }]]]
				set wat_l($j) [linsert $wat_l($j) end [lindex $data [expr { $pair2 + 2 }]]]
				set wat_l($j) [linsert $wat_l($j) end $temp_dis1]
				set wat_l($j) [linsert $wat_l($j) end $temp_dis2]
			}
			#set temp_dis [expr { ($temp_dis1 + $temp_dis2) / 2.0 }]
			set temp_dis $temp_dis1
			set d($j) [expr { $d($j) + $temp_dis }]
		}
	}

	# AVERAGING

	for {set i 1} {$i <= $n_cla} {incr i} {
		set value [expr { $d($i) / $n_frames }]
		set ep1($i) [expr { $ep1($i) / $n_frames }]
		set ep2($i) [expr { $ep2($i) / $n_frames }]
		if { $ep1($i) > $ep2($i) } {
			set min_ep $ep2($i)
		} else { 
			set min_ep $ep1($i)
		}
		puts $g1 " $i $ep1($i) $ep2($i) $min_ep"
		if {$value < 2.5 } {
			puts $g " $i  $value { $wat_l($i) }"
		} 
	}
	close $g
	close $g1
}

proc ester_dis {atom pair1 pair2} {

	# IN THIS PROCEDURE WE WILL CALCULATE THE AVERAGE DISTANCE OF THE METHYL ESTER GROUP OF RING V FROM THE HYDROGENS OF THE WATER BOUND TO THE CENTRAL MAGNESIUM

	set f [open "frame.pdb" "r"]
	set data [read $f]
	close $f

	set k 0
	
	while {[lindex $data $k] != $atom || [lindex $data [expr { $k - 2 }]] != "Mg1"} {
		incr k
	}

	while { [lindex $data $k] != "TER" } {
		if { [lindex $data $k] == "O1D" } {
			set esterx [lindex $data [expr { $k + 3 }]]
			set estery [lindex $data [expr { $k + 4 }]]
			set esterz [lindex $data [expr { $k + 5 }]]
		}
	incr k
	}

	set pair1h1 [expr { [lindex $data $pair1] + 1 }]
	set pair1h2 [expr { [lindex $data $pair1] + 2 }]
	set pair2h1 [expr { [lindex $data $pair2] + 1 }]
	set pair2h2 [expr { [lindex $data $pair2] + 2 }] 

	while { $k < [llength $data] } {
		if { [lindex $data $k] == $pair1h1 && [lindex $data [expr { $k + 3 }]] == [lindex $data [expr { $pair1 + 3 }]]} { 
			set p1h1x [lindex $data [expr { $k + 4 }]]
			set p1h1y [lindex $data [expr { $k + 5 }]]
			set p1h1z [lindex $data [expr { $k + 6 }]]
		}
		if { [lindex $data $k] == $pair1h2 && [lindex $data [expr { $k + 3 }]] == [lindex $data [expr { $pair1 + 3 }]] } { 
			set p1h2x [lindex $data [expr { $k + 4 }]]
			set p1h2y [lindex $data [expr { $k + 5 }]]
			set p1h2z [lindex $data [expr { $k + 6 }]]
		}
		if { [lindex $data $k] == $pair2h1 && [lindex $data [expr { $k + 3 }]] == [lindex $data [expr { $pair2 + 3 }]] } { 
			set p2h1x [lindex $data [expr { $k + 4 }]]
			set p2h1y [lindex $data [expr { $k + 5 }]]
			set p2h1z [lindex $data [expr { $k + 6 }]]
		}
		if { [lindex $data $k] == $pair2h2 && [lindex $data [expr { $k + 3 }]] ==  [lindex $data [expr { $pair2 + 3 }]] } { 
			set p2h2x [lindex $data [expr { $k + 4 }]]
			set p2h2y [lindex $data [expr { $k + 5 }]]
			set p2h2z [lindex $data [expr { $k + 6 }]]
		}
		incr k
	}

	# CALCULATING THE DISTANCE OF PAIR1 AND PAIR2 FROM ESTER GROUP

	set g [open "dummy_ester" "w"]

	# pair1

	set difx [expr { $esterx-$p1h1x }]
	set difx [expr { $difx * $difx }]
	set dify [expr { $estery-$p1h1y }]
	set dify [expr { $dify * $dify }]
	set difz [expr { $esterz-$p1h1z }]
	set difz [expr { $difz * $difz }]

	set dis1 [expr { sqrt($difx + $dify + $difz) }]

	set difx [expr { $esterx-$p1h2x }]
	set difx [expr { $difx * $difx }]
	set dify [expr { $estery-$p1h2y }]
	set dify [expr { $dify * $dify }]
	set difz [expr { $esterz-$p1h2z }]
	set difz [expr { $difz * $difz }]

	set dis2 [expr { sqrt($difx + $dify + $difz) }]

	if {$dis1 < $dis2} {
		puts $g "{ $dis1 }"
	} else { 
		puts $g "{ $dis2 }"
	}

	# pair2

	set difx [expr { $esterx-$p2h1x }]
	set difx [expr { $difx * $difx }]
	set dify [expr { $estery-$p2h1y }]
	set dify [expr { $dify * $dify }]
	set difz [expr { $esterz-$p2h1z }]
	set difz [expr { $difz * $difz }]

	set dis1 [expr { sqrt($difx + $dify + $difz) }]

	set difx [expr { $esterx-$p2h2x }]
	set difx [expr { $difx * $difx }]
	set dify [expr { $estery-$p2h2y }]
	set dify [expr { $dify * $dify }]
	set difz [expr { $esterz-$p2h2z }]
	set difz [expr { $difz * $difz }]

	set dis2 [expr { sqrt($difx + $dify + $difz) }]

	if {$dis1 < $dis2} {
		puts $g "{ $dis1 }"
	} else { 
		puts $g "{ $dis2 }"
	}
	close $g
}
	
proc mg_pair {} {
	
	# IN THIS PROCEDURE WE WILL CALCULATE THE NEAREST NEIGBOUR OF EACH MG AND CALCULATE THE VORNOI POLYHEDRA OF EACH ATOM IN THAT RESIDUE WRT TO CENTRAL MG OF THE RESIDUE UNDER CONSIDERATION

	package require math::linearalgebra 

	set n_cla [input2]

	set num [open "[input4]" "r"]	
	set natom [read $num]
	close $num

	set k 0
	while { [lindex $natom $k] != "%FORMAT(10I8)" } {
		incr k
	}
	incr k
	set n_atoms [lindex $natom $k]
	set num_at_pres [input3]

	for {set i 0} {$i <= $n_cla} {incr i} {
		set d($i) 0.0
		set l($i) ""
		set n_dl($i) ""
		for {set j 0} {$j <= $n_cla} {incr j} {
			set dpair($i,$j) 0.0
		}
	}
	
	set m [open "mg_mg_distance" "w"]

	set mp [open "mg_distance_all" "w"]

	set mgc [open "mg_coord" "w"]

	set start_frame [input5]
	set end_frame [input6]
	set step [input7]
	set n_frames [expr { (($end_frame - $start_frame)/$step) + 1 }]

	for {set fr $start_frame} {$fr <= $end_frame} {incr fr $step} {

		# EXECUTING CPPTRAJ

		set inp [open "ai.dat" "w"]
		puts $inp "trajin [cpptraj_traj] $fr $fr"
		puts $inp "trajout frame.mdcrd mdcrd nobox"
		puts $inp "trajout frame.rst rst"
		puts $inp "go"
		close $inp

		set inp [open "ai.dat" "r"]

		exec cpptraj -p [input4] -i ai.dat 

		close $inp

		set f1 [open "frame.rst" "r"]
		set data1 [read $f1]
		close $f1

		set box_x [lindex $data1 [expr { [llength $data1] - 3 }]]
		set box_y [lindex $data1 [expr { [llength $data1] - 2 }]]
		set box_z [lindex $data1 [expr { [llength $data1] - 1 }]]

		puts "				**** FRAME :: $fr"
		set f [open "$fr" "w"]
		set coord [com]
		puts $f "{$coord}"
		close $f

		# CALCULATING THE DISTANCE OF EACH RESIDUE FROM OTHER RESIDUES

		set f [open "$fr" "r"]
		set data [read $f]
		close $f
		for {set i 1} {$i <= $n_cla} {incr i} {
			set res1x [lindex $data 0 [expr { $i - 1 }] 0]
			set res1y [lindex $data 0 [expr { $i - 1 }] 1]
			set res1z [lindex $data 0 [expr { $i - 1 }] 2]
		
			puts $mgc "$res1x $res1y $res1z"

			set dis_temp 1000.0

			for {set j 1} {$j <= $n_cla} {incr j} {
				if {$j != $i } {
					set res2x [lindex $data 0 [expr { $j - 1 }] 0]
					set res2y [lindex $data 0 [expr { $j - 1 }] 1]
					set res2z [lindex $data 0 [expr { $j - 1 }] 2]

					# CALCULATING THE DISTANCE BETWEEN THE TWO RESIDUES

					set difx [expr { abs($res1x-$res2x) }]
					if { $difx > [expr { $box_x / 2.0 }] } {
						set difx [expr { $difx - $box_x }]
					}
					set difx [expr { $difx * $difx }]

					set dify [expr { abs($res1y-$res2y) }]
					if { $dify > [expr { $box_y / 2.0 }] } {
						set dify [expr { $dify - $box_y }]
					}
					set dify [expr { $dify * $dify }]

					set difz [expr { abs($res1z-$res2z) }]
					if { $difz > [expr { $box_z / 2.0 }] } {
						set difz [expr { $difz - $box_z }]
					}
					set difz [expr { $difz * $difz }]

					set dis [expr { sqrt($difx + $dify + $difz) }]

					set dpair($i,$j) [expr { $dpair($i,$j) + $dis }]

					if { $dis < $dis_temp } {
						set dis_temp $dis
						set pair $j
					}			
				}
			}
			set d($i) [expr { $d($i) + $dis_temp }]	
			if { [expr { $fr % 10 }] == 0 } {
				set l($i) [linsert $l($i) end $pair]
				set n_dl($i) [linsert $n_dl($i) end $dis_temp]
			}
		}
	} 
	close $mgc
	# AVERAGING

	for {set i 1} {$i <= $n_cla} {incr i} {
		set value [expr { $d($i) / $n_frames }]
		puts $m " $i $value { $l($i) $n_dl($i) }"
		for {set j 1} {$j <= $n_cla} {incr j} {
			if { $j != $i } {
				set value [expr {$dpair($i,$j) / $n_frames }]
				puts $mp " $i $j $value "
			}
		}
	}
	close $m
}
	
proc ester_mg_dis {} {

	# DETERMINING THE NEAREST WATER NEIGBOUR TO THE MG IN EACH RESIDUES

	set n_cla [input2]	

	for {set i 1} {$i <= $n_cla} {incr i} {
		set d($i) 0.0
		set wat_l($i) ""
	}
	
	set g [open "ESTER_MG_dis" "w"]

	set start_frame [input5]
	set end_frame [input6]
	set step [input7]
	set n_frames [expr { (($end_frame - $start_frame)/$step) + 1 }]

	for {set fr $start_frame} {$fr <= $end_frame} {incr fr $step} {

		# EXECUTING CPPTRAJ

		set inp [open "ai.dat" "w"]
		puts $inp "trajin [cpptraj_traj] $fr $fr"
		puts $inp "trajout frame.pdb pdb nobox"
		puts $inp "trajout frame.rst rst"
		puts $inp "go"
		close $inp

		set inp [open "ai.dat" "r"]

		exec cpptraj -p [input4] -i ai.dat 

		close $inp

		set f [open "frame.pdb" "r"]
		set data [read $f]
		close $f
	
		set f1 [open "frame.rst" "r"]
		set data1 [read $f1]
		close $f1

		set box_x [lindex $data1 [expr { [llength $data1] - 3 }]]
		set box_y [lindex $data1 [expr { [llength $data1] - 2 }]]
		set box_z [lindex $data1 [expr { [llength $data1] - 1 }]]

		for {set j 1} {$j <= $n_cla} {incr j} {

			set d_temp 0.0

			puts "		****** FRAME $fr :: step $j of $n_cla :: of [input6] frames   ******		"

			set k 0
			while { [lindex $data $k] != $j || [lindex $data [expr { $k - 2 }]] != "O1A" } {
				incr k
			}
			set O1Dx [lindex $data [expr { $k + 1 }]]
			set O1Dy [lindex $data [expr { $k + 2 }]]
			set O1Dz [lindex $data [expr { $k + 3 }]]

			set k 0

			set temp_dis1 1000.0
			set temp_dis2 1000.0
			set pair1 1000000

			while { [lindex $data $k] != "END" } {
				if { [lindex $data $k] == "Mg1" && [lindex $data [expr { $k + 1 }]] == "MOL" } {
					set ox [lindex $data [expr { $k + 3 }]]
					set oy [lindex $data [expr { $k + 4 }]]
					set oz [lindex $data [expr { $k + 5 }]]

					set difx [expr { abs($O1Dx-$ox) }]
					if { $difx > [expr { $box_x / 2.0 }] } {
						set difx [expr { $difx - $box_x }]
					}
					set difx [expr { $difx * $difx }]

					set dify [expr { abs($O1Dy-$oy) }]
					if { $dify > [expr { $box_y / 2.0 }] } {
						set dify [expr { $dify - $box_y }]
					}
					set dify [expr { $dify * $dify }]
		
					set difz [expr { abs($O1Dz-$oz) }]
					if { $difz > [expr { $box_z / 2.0 }] } {
						set difz [expr { $difz - $box_z }]
					}
					set difz [expr { $difz * $difz }]

					set dis [expr { sqrt($difx + $dify + $difz) }]
		
					if { $dis < $temp_dis1 } {
							set temp_dis2 $temp_dis1
							set temp_dis1 $dis
							set pair2 $pair1
							set pair1 $k					
					} elseif { $dis < $temp_dis2 } {
							set temp_dis2 $dis
							set pair2 $k
					} 
				}
				incr k 
			}
			if { [expr { $fr % 10 }] == 0 } {
				set wat_l($j) [linsert $wat_l($j) end [lindex $data [expr { $pair1 + 2 }]]]
				set wat_l($j) [linsert $wat_l($j) end [lindex $data [expr { $pair2 + 2 }]]]
				set wat_l($j) [linsert $wat_l($j) end $temp_dis1]
				set wat_l($j) [linsert $wat_l($j) end $temp_dis2]
			}
			#set temp_dis [expr { ($temp_dis1 + $temp_dis2) / 2.0 }]
			set temp_dis $temp_dis1
			set d($j) [expr { $d($j) + $temp_dis }]
		}
	}

	# AVERAGING

	for {set i 1} {$i <= $n_cla} {incr i} {
		set value [expr { $d($i) / $n_frames }]
		puts $g " $i  $value { $wat_l($i) }" 
	}
	close $g
}

proc delete_files {} {
	
	# THIS PROCEDURE WILL DELETE ALL THE FILES REQUIRED DURING CALCULATIONS ONLY

	for {set i 0} {$i < [input2]} {incr i} {
		file delete $i
		file delete $i.avg
	}

	set start_frame [input5]
	set end_frame [input6]
	for {set fr $start_frame} {$fr <= $end_frame} {incr fr} {
		file delete $fr
	}

	file delete dummy
	file delete dummy_ester
	#file delete frame.pdb 
	file delete frame.mdcrd
	file delete ai.dat

}

proc input1 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 0]
}

proc input2 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 1]
}

proc input3 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 2]
}

proc input4 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 3]
}

proc input5 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 5 0]
}

proc input6 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 5 1]
}

proc input7 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 5 2]
} 

proc cpptraj_traj {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 6]
}

proc output1 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 4 0]
}

proc output2 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 4 1]
}
#plane_N_static
#plane_N
#water_distance
#ester_mg_dis
mg_pair
delete_files























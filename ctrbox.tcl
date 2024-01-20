# centres the waterbox regarding the geometry

proc ctrbox { in_psf in_pdb out_pfx } {

	resetpsf
	readpsf $in_psf
	coordpdb $in_pdb

	mol load psf $in_psf pdb $in_pdb

	set all [ atomselect top all ]
	set minmax [ measure minmax $all ]

	foreach {min max} $minmax { break }
	foreach {xmin ymin zmin} $min { break }
	foreach {xmax ymax zmax} $max { break }

	set mx [ expr $xmin + abs(($xmin - $xmax) / 2) ]
	set my [ expr $ymin + abs(($ymin - $ymax) / 2) ]
	set mz [ expr $zmin + abs(($zmin - $zmax) / 2) ]

	$all moveby [ vecsub {0 0 0} [list $mx $my $mz] ]
	
	foreach atom [$all get {segid resid name x y z}] {
		foreach {segid resid name x y z} $atom { break }
		coord $segid $resid $name [list $x $y $z]
	}
	
	writepsf $out_pfx.psf 
	writepdb $out_pfx.pdb

	set minmax [ measure minmax $all ]
	foreach {min max} $minmax { break }
	foreach {xmin ymin zmin} $min { break }

	  set vec1 [ expr 2 * abs($xmin) + 0.5 ]
	  set vec2 [ expr 2 * abs($ymin) + 0.5 ]
	  set vec3 [ expr 2 * abs($zmin) + 0.5 ]
        
        set fp [ open "$out_pfx.pdb" a+ ]
          puts $fp "REMARK cellBasisVector1	$vec1	0.0	0.0"
	  puts $fp "REMARK cellBasisVector2 	0.0	$vec2	0.0"
	  puts $fp "REMARK cellBasisVector3 	0.0	0.0	$vec3"
	  puts $fp "REMARK cellOrigin		0.0	0.0	0.0"

        close $fp

	mol delete top
}

package require pbctools
# Load trajectory
mol new "../leut-f108y-apo1_full.psf" type psf waitfor all
animate read dcd  "/gs-scratch/csb/bradshaw/3TT1_F108Y/Run_1/MD/run_analysis/2us_full_run1.dcd" waitfor all
pbc wrap -all -compound res -center com -centersel protein

# Distances chosen for 'opening' of vestibule as per Zomot E. et al., DOI: 10.1074/jbc.M114.617555
# TM1b-TM10 distance
set out [open "TM1b-TM10_distance.dat" a];
puts $out "# Distance V33(C) - D401 (Ca)"

set idxs [[atomselect top "protein and resid 33 and name C"] list]
lappend idxs [[atomselect top "protein and resid 401 and name CA"] list]

set num_steps [molinfo top get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    set currdist [measure bond $idxs frame $frame]
    puts $out [format "%.5f" $currdist]
    }
close $out

# TM6a-TM10 distance
set out2 [open "TM6a-TM10_distance.dat" a];
puts $out2 "# Distance I245(C) - I410(Ca)"

set idxs [[atomselect top "protein and resid 245 and name C"] list]
lappend idxs [[atomselect top "protein and resid 410 and name CA"] list]

set num_steps [molinfo top get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    set currdist [measure bond $idxs frame $frame]
    puts $out2 [format "%.5f" $currdist]
    }
close $out2

# Residue Ca distances to measure occlusion, as per Toni
# Dists = Y108-L25, Y108-F253, D404-R30
set out3 [open "Ca-Ca_distances_measuring_occlusion.dat" a];
puts $out3 "#Y108-L25 Y108-F253 D404-R30"

set tyr [[atomselect top "protein and resid 108 and name CA"] list]
set leu [[atomselect top "protein and resid 25 and name CA"] list]
set phe [[atomselect top "protein and resid 253 and name CA"] list]
set asp [[atomselect top "protein and resid 404 and name CA"] list]
set arg [[atomselect top "protein and resid 30 and name CA"] list]

set num_steps [molinfo top get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    set dists [format "%.5f" [measure bond "$tyr $leu" frame $frame]]
    lappend dists [format "%.5f" [measure bond "$tyr $phe" frame $frame]]
    lappend dists [format "%.5f" [measure bond "$asp $arg" frame $frame]]
    puts $out3 $dists
}
close $out3

exit

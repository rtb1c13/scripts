package require pbctools
# Load trajectory
mol new "../leut-f108y-apo1_full.psf" type psf waitfor all
animate read dcd  "/gs-scratch/csb/bradshaw/3TT1_F108Y/Run_1/MD/run_analysis/2us_full_run1.dcd" waitfor all
pbc wrap -all -compound res -center com -centersel protein

#Na+ 1 coordinated to N27sc,T254sc,N286sc,A22bb,T254bb
set out [open "Na+1_distances.dat" a];
puts $out "#N27sc  T254sc  N286sc  A22bb   T254bb"

set sod [[atomselect top "segid SODI and resid 1"] list]
set idxs [[atomselect top "protein and resid 27 and name OD1"] list]
lappend idxs [[atomselect top "protein and resid 254 and name OG1"] list]
lappend idxs [[atomselect top "protein and resid 286 and name OD1"] list]
lappend idxs [[atomselect top "protein and resid 22 and name O"] list]
lappend idxs [[atomselect top "protein and resid 254 and name O"] list]


set num_steps [molinfo top get numframes]
set dists ""
for {set frame 0} {$frame < $num_steps} {incr frame} {
    set dists ""
    foreach a $idxs {
        set currdist [measure bond "$sod $a" frame $frame]
        lappend dists [format "%.5f" $currdist]
    }
    puts $out $dists
    unset dists
}
close $out

#Na+ 2 coordinated to T354sc,S355sc,G20bb,V23bb,A351bb
set out2 [open "Na+2_distances.dat" a];
puts $out2 "#T354sc S355sc  G20sc   V23bb   A351bb"

set sod [[atomselect top "segid SODI and resid 2"] list]
set idxs [[atomselect top "protein and resid 354 and name OG1"] list]
lappend idxs [[atomselect top "protein and resid 355 and name OG"] list]
lappend idxs [[atomselect top "protein and resid 20 and name O"] list]
lappend idxs [[atomselect top "protein and resid 23 and name O"] list]
lappend idxs [[atomselect top "protein and resid 351 and name O"] list]


set num_steps [molinfo top get numframes]
set dists ""
for {set frame 0} {$frame < $num_steps} {incr frame} {
    set dists ""
    foreach a $idxs {
        set currdist [measure bond "$sod $a" frame $frame]
        lappend dists [format "%.5f" $currdist]
    }
    puts $out $dists
    unset dists
}
close $out2

exit

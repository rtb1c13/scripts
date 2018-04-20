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

exit

# Adapted from $ANU/from_pacific/leut-f108y-2/sasa_tm.tcl
#set j [lindex $argv 0]

mol new ../../leut-f108y-apo1_full.psf type psf

set skip 1
#1 to 16 for segs 1-3
for {set res 1} {$res < 17} {incr res 1} {
puts $res
set outfile1 [open "res$res.dat" a];
set prtlip [atomselect top "( segid PROT or resname DMPC)"]
set protein [atomselect top "( segid PROT and resid $res )"]
set bb [atomselect top "( segid PROT and resid $res and ( name CA or name C or name N or name O ))"]
set sc [atomselect top "( segid PROT and resid $res and not ( name CA or name C or name N or name O))"]
# sasa calculation loop
animate read dcd ../2us_full_run1.dcd  waitfor all
set num_steps [molinfo top get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    molinfo top set frame $frame
    set sasa [measure sasa 1.4 $prtlip -restrict $protein]
    set sasa1 [measure sasa 1.4 $prtlip -restrict $bb]
    set sasa2 [measure sasa 1.4 $prtlip -restrict $sc]
    puts $outfile1 "$sasa $sasa1 $sasa2"
}
close $outfile1
}


exit



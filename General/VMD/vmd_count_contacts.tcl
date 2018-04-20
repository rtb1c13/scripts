package require pbctools
# Load trajectory
mol new "../leut-f108y-apo1_full.psf" type psf waitfor all
animate read dcd  "2us_full_run1.dcd" waitfor all
pbc wrap -all -compound res -center com -centersel protein

# Output
set out [open "DMPC_within_6.5A.dat" w]

# Count
set sel [atomselect top "same residue as (resname DMPC and within 6.5 of protein)"]
set num_steps [molinfo top get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
  $sel frame $frame
  $sel update
  set n [$sel num]
  # 118 = atoms in DMPC
  puts $out "[expr $n / 118]"
}
close $out
quit

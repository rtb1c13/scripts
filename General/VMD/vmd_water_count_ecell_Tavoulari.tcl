# Water count in LeuT extracellular vestibule volume defined as:
# Lipid phosphate COM = 0 in z-axis
# Protein COM = 0 in x/y directions
# Waters within box of dimensions (20, 15, 16)A, centered at (0, -2.5, 2)
# As per: Tavoulari S. et al., DOI: 10.1074/jbc.M115.692012

package require pbctools

# Key bindings for rotation 
user add key q {rotate x by 90}
user add key w {rotate y by 90}
user add key e {rotate z by 90}

# Load trajectory
mol new "../leut-f108y-apo1_full.psf" type psf waitfor all
animate read dcd  "/media/bradshaw/RTB_externa/lara/WORK/LeuT/Anu_files/3TT1_f108y/Run_1/MD/run_analysis/2us_full_run1.dcd" beg 15000 end -1 waitfor all
pbc wrap -all -compound res -center com -centersel protein


set outfile1 [ open "watcount_ecell_Tavoulari_run1.dat" a ];
set num_steps [ molinfo top get numframes ] 

#set K to whatever last frame number you have
for {set k 0 } {$k < $num_steps} { incr k } {

#    animate read dcd $dcd  beg $k end $k
    # Calc current offsets of COM
    set sel [ atomselect top "resname DMPC and name P" frame $k ]
    set com [ measure center $sel weight mass ]
    set z1 [expr 0 - [ lindex $com 2 ]]
    set sel [ atomselect top "protein" frame $k] 
    set com [ measure center $sel weight mass ]
    set x2 [expr 0 - [ lindex $com 0 ]]
    set y2 [expr 0 - [ lindex $com 1 ]]
    set sel [ atomselect top "all" frame $k]
    # Translate origin to COM of protein in x/y, membrane phosphates in z
    $sel move [ transoffset "$x2 $y2 $z1"] 

    set minx -10
    set maxx 10
    set miny -10
    set maxy 5
    set minz -6 
    set maxz 10

    draw materials off
    draw color blue
    draw line "$minx $miny $minz" "$maxx $miny $minz" width 3
    draw line "$minx $miny $minz" "$minx $maxy $minz" width 3
    draw line "$minx $miny $minz" "$minx $miny $maxz" width 3

    draw line "$maxx $miny $minz" "$maxx $maxy $minz" width 3
    draw line "$maxx $miny $minz" "$maxx $miny $maxz" width 3

    draw line "$minx $maxy $minz" "$maxx $maxy $minz" width 3
    draw line "$minx $maxy $minz" "$minx $maxy $maxz" width 3

    draw line "$minx $miny $maxz" "$maxx $miny $maxz" width 3
    draw line "$minx $miny $maxz" "$minx $maxy $maxz" width 3

    draw line "$maxx $maxy $maxz" "$maxx $maxy $minz" width 3
    draw line "$maxx $maxy $maxz" "$minx $maxy $maxz" width 3
    draw line "$maxx $maxy $maxz" "$maxx $miny $maxz" width 3
    
    set selwat [atomselect top "same residue as (resname TIP3 and ( x > $minx and x < $maxx) and ( y > $miny and y < $maxy) and ( z > $minz and z < $maxz))" frame $k ]
 
#WATER RESIDUES
    #set watres [llength [lsort -integer -unique [$selwat get residue]]]
    #puts $outfile1 "$watres"

    set n [$selwat num]
    # 3 = atoms in TIP3P
    puts $outfile1 "[expr $n / 3]"
}

close $outfile1
exit

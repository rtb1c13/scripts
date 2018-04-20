# Water count in LeuT extracellular vestibule volume defined as:
# Within 26A of F259-CB
# Not within 5A of lipid
# Z coordinate j between F259-CB-z < j < F259-CB-z + 23
# As per: C. Zhao et al., Biophys J, 2012, 103 (5), 878-888

package require pbctools
#set xray "/u/bradshaw/WORK/LeuT/Anu_files/3TT1_f108y/Xtal/3tt1_molprobity_mod.pdb"; # Xtal structure, truncated from resid 1-4, 130-134, after molprobity
#mol new ${xray}
#mol delrep 0 top
# Load trajectory
mol new "../leut-f108y-apo1_full.psf" type psf waitfor all
animate read dcd  "/media/bradshaw/RTB_externa/lara/WORK/LeuT/Anu_files/3TT1_f108y/Run_1/MD/run_analysis/2us_full_run1.dcd" beg 15000 end -1 waitfor all
pbc wrap -all -compound res -center com -centersel protein

set outfile1 [ open "watcount_run1.dat" a ];
set num_steps [ molinfo top get numframes ] 

# Align hash domain to first frame
set reference [ atomselect top "protein and backbone and (resid 88 to 124 or resid 166 to 185 or resid 337 to 371 or resid 375 to 393)" frame 0 ]  
set compare [ atomselect top "protein and backbone and (resid 88 to 124 or resid 166 to 185 or resid 337 to 371 or resid 375 to 393)" ]  

#set K to whatever last frame number you have
for {set k 0 } {$k < $num_steps} { incr k } {

#    animate read dcd $dcd  beg $k end $k
    $compare frame $k
    set trans_mat [ measure fit $compare $reference ]
    [ atomselect top "all" frame $k ] move $trans_mat

    # Get Z-coordinate of 259-CB
    set selcb [ atomselect top "protein and resid 259 and name CB" frame $k ]
    set minz [ $selcb get {z} ]
    set maxz [ expr $minz + 23 ]

    set selwat [atomselect top "{same residue as ((resname TIP3 and (not within 5 of resname DMPC)) and ( z > $minz and z < $maxz) and ( within 26 of ( protein and resid 259 and name CB )))}" frame $k ]
 
#WATER RESIDUES
    #set watres [llength [lsort -integer -unique [$selwat get residue]]]
    #puts $outfile1 "$watres"

    set n [$selwat num]
    # 3 = atoms in TIP3P
    puts $outfile1 "[expr $n / 3]"
}

close $outfile1
exit

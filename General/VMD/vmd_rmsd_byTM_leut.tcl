# RMSD of individual LeuT helices printed to file 

### Specify trajectory & X-ray for analysis ###
set traj_dir "/media/bradshaw/RTB_externa/lara/WORK/LeuT/Anu_files/3TT1_f108y/Run_1/MD/"
set xpsf_file "../3TT1_protonly.psf"; 
set xray "/u/bradshaw/WORK/LeuT/Anu_files/3TT1_f108y/Run_1/MD/rsts/leut-f108y-apo1_full.emin.crd"; # Xtal structure, after setup, griffin, minimisation, before MD
###

# Load X-ray (Charmm coords)
mol new ${xray} type cor
mol rename top {X-ray}

# Load trajectory
mol new ${xpsf_file} type psf waitfor all
#foreach fname [glob ${traj_dir}*.dcd] {
animate read dcd  "${traj_dir}run_analysis/2us_protonly_run1.dcd" waitfor all
#}

# RMSD procedure for all frames
proc rmsd_allframes {{mol top} refid refsel molsel} {
    set reference [atomselect $refid $refsel frame 0]
    set compare [atomselect $mol $molsel] 
    set all_ref [atomselect $refid "all" frame 0]
    set all_compare [atomselect $mol "all"]
    set num_steps [molinfo top get numframes]
    for {set frame 0} {$frame < $num_steps} {incr frame} {
        $compare frame $frame
        $all_compare frame $frame
        set trans_mat [measure fit $compare $reference]
        $all_compare move $trans_mat
        set rms($frame) [measure rmsd $compare $reference]
    }
    return [array get rms]
}
# RMSD procedure for one frame
proc rmsd_currframe {{mol top} refid refsel molsel currframe} {
    set reference [atomselect $refid $refsel frame 0]
    set compare [atomselect $mol $molsel] 
    set all_ref [atomselect $refid "all" frame 0]
    set all_compare [atomselect $mol "all"]
    $compare frame $currframe
    $all_compare frame $currframe
    set trans_mat [measure fit $compare $reference]
    $all_compare move $trans_mat
    set rms [measure rmsd $compare $reference]
    return $rms
}

# RMSD procedure for one frame without fitting
proc rmsd_nofit {{mol top} refid refsel molsel currframe} {
    set reference [atomselect $refid $refsel frame 0]
    set compare [atomselect $mol $molsel] 
    set all_ref [atomselect $refid "all" frame 0]
    set all_compare [atomselect $mol "all"]
    $compare frame $currframe
    $all_compare frame $currframe
    set rms [measure rmsd $compare $reference]
    return $rms
}

# Alignments & RMSDs
set num_steps [molinfo top get numframes]
set outfile1 [ open "Byhelix_RMSD_fitted_to_self.dat" a ];
puts $outfile1 "#TMs defined identically to Adhikary et al., DOI: 10.1073/pnas.1613293114"
puts $outfile1 "#TM1    TM2     TM3     TM4     TM5     TM6     TM7     TM8     TM9     TM10    TM11    TM12    IL1     EL2     EL3     EL4"
for {set frame 0} {$frame < $num_steps} {incr frame} {
    # TM1
    set currrms [rmsd_currframe top 0 "backbone and (resid 12 to 38)" "segid PROT and backbone and (resid 12 to 38)" $frame]
    set allrms [format "%7.5f" $currrms] 
    # TM2
    set currrms [rmsd_currframe top 0 "backbone and (resid 42 to 72)" "segid PROT and backbone and (resid 42 to 72)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM3
    set currrms [rmsd_currframe top 0 "backbone and (resid 88 to 124)" "segid PROT and backbone and (resid 88 to 124)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM4
    set currrms [rmsd_currframe top 0 "backbone and (resid 166 to 184)" "segid PROT and backbone and (resid 166 to 184)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM5
    set currrms [rmsd_currframe top 0 "backbone and (resid 192 to 215)" "segid PROT and backbone and (resid 192 to 215)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM6
    set currrms [rmsd_currframe top 0 "backbone and (resid 241 to 268)" "segid PROT and backbone and (resid 241 to 268)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM7
    set currrms [rmsd_currframe top 0 "backbone and (resid 276 to 306)" "segid PROT and backbone and (resid 276 to 306)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM8
    set currrms [rmsd_currframe top 0 "backbone and (resid 337 to 369)" "segid PROT and backbone and (resid 337 to 369)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM9
    set currrms [rmsd_currframe top 0 "backbone and (resid 375 to 395)" "segid PROT and backbone and (resid 375 to 395)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM10
    set currrms [rmsd_currframe top 0 "backbone and (resid 399 to 427)" "segid PROT and backbone and (resid 399 to 427)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM11
    set currrms [rmsd_currframe top 0 "backbone and (resid 446 to 476)" "segid PROT and backbone and (resid 446 to 476)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM12
    set currrms [rmsd_currframe top 0 "backbone and (resid 483 to 510)" "segid PROT and backbone and (resid 483 to 510)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # IL1
    set currrms [rmsd_currframe top 0 "backbone and (resid 77 to 84)" "segid PROT and backbone and (resid 77 to 84)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # EL2
    set currrms [rmsd_currframe top 0 "backbone and (resid 137 to 153)" "segid PROT and backbone and (resid 137 to 153)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # EL3
    set currrms [rmsd_currframe top 0 "backbone and (resid 223 to 231)" "segid PROT and backbone and (resid 223 to 231)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # EL4
    set currrms [rmsd_currframe top 0 "backbone and (resid 308 to 317 or resid 320 to 331)" "segid PROT and backbone and (resid 308 to 317 or resid 320 to 331)" $frame]
    lappend allrms [format "%7.5f" $currrms]

    puts $outfile1 $allrms
}
close $outfile1

# RMSD after fitting to hash region
set outfile2 [ open "Byhelix_RMSD_fitted_to_hash.dat" a ];
puts $outfile2 "#TMs defined identically to Adhikary et al., DOI: 10.1073/pnas.1613293114"
puts $outfile2 "#TM1    TM2     TM3     TM4     TM5     TM6     TM7     TM8     TM9     TM10    TM11    TM12    IL1     EL2     EL3     EL4"
# Fit
rmsd_allframes top 0 "backbone and (resid 88 to 124 or resid 166 to 185 or resid 337 to 371 or resid 375 to 393)" "segid PROT and backbone and (resid 88 to 124 or resid 166 to 185 or resid 337 to 371 or resid 375 to 393)"
for {set frame 0} {$frame < $num_steps} {incr frame} {
    # TM1
    set currrms [rmsd_nofit top 0 "backbone and (resid 12 to 38)" "segid PROT and backbone and (resid 12 to 38)" $frame]
    set allrms [format "%7.5f" $currrms] 
    # TM2
    set currrms [rmsd_nofit top 0 "backbone and (resid 42 to 72)" "segid PROT and backbone and (resid 42 to 72)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM3
    set currrms [rmsd_nofit top 0 "backbone and (resid 88 to 124)" "segid PROT and backbone and (resid 88 to 124)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM4
    set currrms [rmsd_nofit top 0 "backbone and (resid 166 to 184)" "segid PROT and backbone and (resid 166 to 184)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM5
    set currrms [rmsd_nofit top 0 "backbone and (resid 192 to 215)" "segid PROT and backbone and (resid 192 to 215)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM6
    set currrms [rmsd_nofit top 0 "backbone and (resid 241 to 268)" "segid PROT and backbone and (resid 241 to 268)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM7
    set currrms [rmsd_nofit top 0 "backbone and (resid 276 to 306)" "segid PROT and backbone and (resid 276 to 306)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM8
    set currrms [rmsd_nofit top 0 "backbone and (resid 337 to 369)" "segid PROT and backbone and (resid 337 to 369)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM9
    set currrms [rmsd_nofit top 0 "backbone and (resid 375 to 395)" "segid PROT and backbone and (resid 375 to 395)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM10
    set currrms [rmsd_nofit top 0 "backbone and (resid 399 to 427)" "segid PROT and backbone and (resid 399 to 427)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM11
    set currrms [rmsd_nofit top 0 "backbone and (resid 446 to 476)" "segid PROT and backbone and (resid 446 to 476)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # TM12
    set currrms [rmsd_nofit top 0 "backbone and (resid 483 to 510)" "segid PROT and backbone and (resid 483 to 510)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # IL1
    set currrms [rmsd_nofit top 0 "backbone and (resid 77 to 84)" "segid PROT and backbone and (resid 77 to 84)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # EL2
    set currrms [rmsd_nofit top 0 "backbone and (resid 137 to 153)" "segid PROT and backbone and (resid 137 to 153)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # EL3
    set currrms [rmsd_nofit top 0 "backbone and (resid 223 to 231)" "segid PROT and backbone and (resid 223 to 231)" $frame]
    lappend allrms [format "%7.5f" $currrms]
    # EL4
    set currrms [rmsd_nofit top 0 "backbone and (resid 308 to 317 or resid 320 to 331)" "segid PROT and backbone and (resid 308 to 317 or resid 320 to 331)" $frame]
    lappend allrms [format "%7.5f" $currrms]

    puts $outfile2 $allrms
}
close $outfile2
exit


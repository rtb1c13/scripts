# Source into VMD. Adapted from VMD mailing list.
# Use as: draw arrow "index 1000" { lenx leny lenz }


# ColorID 4 = yellow (highly visible)
graphics top color 4

proc vmd_draw_arrow {mol selstring veclens} {
     # an arrow is made of a cylinder and a cone
     set sel [atomselect top $selstring]
     if {[$sel num] != 1} {
         error "draw_arrow: '$selstring' must select 1 atom"
     }
     # Get coords
     lassign [$sel get {x y z}] start
     set end [vecadd $start $veclens]
     set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
     graphics $mol cylinder $start $middle radius 0.10
     graphics $mol cone $middle $end radius 0.25
}


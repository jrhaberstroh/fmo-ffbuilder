draw delete all
draw color red
animate goto 0
display resetview
set nm_index [lindex $argv 0]
for {set i 0} {$i < 140} {incr i} {
    set sel_coord [atomselect top "index $i" frame 0]
    set coord [lindex [$sel_coord get {x y z}] 0]
    set sel_ev [atomselect top "index $i" frame $nm_index]
    set eigvec [lindex [$sel_ev get {x y z}] 0]

    set cyl_scale .4
    set cone_scale [expr $cyl_scale * 1.4]
    set cyl_end [vecadd $coord [vecscale $eigvec $cyl_scale]]
    set cone_end [vecadd $coord [vecscale $eigvec $cone_scale]]

    draw cylinder $coord $cyl_end radius 0.05
    draw cone $cyl_end $cone_end radius 0.15
}

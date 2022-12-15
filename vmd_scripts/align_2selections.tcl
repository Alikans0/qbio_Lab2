#example: rmsd 0 1 "name BB and resid 346 to 355 and segid PROD" "name BB and resid 344 to 353" "name BB and resid 5 to 213 223 to 301 and segid PROG PROH" "name BB and resid 4 to 212 222 to 300 and index < 683" "test.dat"

#molid1 = static frame (Reference)
#molid2 = trajectory
proc align2sels {molid1 molid2 frameRef seleAlign1 seleAlign2} {
    set refAlign1 [atomselect $molid1 $seleAlign1 frame $frameRef]
    set num_steps [molinfo $molid2 get numframes]
    set actuelAlign2 [atomselect $molid2 $seleAlign2 frame 0]
    set all [atomselect $molid2 all]
    for {set i 0} {$i < $num_steps} {incr i} {
        $actuelAlign2 frame $i
        $actuelAlign2 update
        $all frame $i
        $all update
        set trans_mat [measure fit $actuelAlign2 $refAlign1]
    	puts "$trans_mat"
        $all move $trans_mat
    }
}


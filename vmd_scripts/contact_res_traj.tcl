set ldlaResidList [list 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65]
mol new "/home/ali/RELAXIN/CG_MD/full_traj1/reimaged.gro" type gro waitfor all
mol addfile "/home/ali/RELAXIN/CG_MD/full_traj1/clusters_simu.trr" type trr waitfor all
#set ldlaResidList [list 20]
foreach ldlaResid $ldlaResidList {
    set fo [open contact_$ldlaResid.dat w]
    set res "resid $ldlaResid"
    set sel [atomselect top "name BB"]
    set residueList [$sel get resid]
    set countList {}
    set resIdList {}
    foreach resId $residueList {
        lappend resIdList [[atomselect top "resid $resId and name BB"] get resid]
        lappend countList 0
    }
    #puts $countList
    #puts $resIdList
    set numFrame [molinfo top get numframes]
    for {set frame 0} {$frame < $numFrame} {incr frame} {
        set contact [atomselect top "(within 8.0 of $res) and resid > 100 and name BB" frame $frame]
        $contact frame $frame
        $contact update
        set contactList [$contact get resid]
        foreach contactRes $contactList {
            set x [expr [lindex [lindex $countList $contactRes] 0] + 1]
            set countList [lreplace $countList $contactRes $contactRes $x]
        }
    }
    set listLength [llength $countList]
    for {set i 0} {$i < $listLength} {incr i} {
        set count [lindex $countList $i]
        if {$count != 0} {
            set count [expr ($count / 10.0)*100.0]
            set resid [lindex $resIdList $i]
            set resname [[atomselect top "resid $resid and name BB"] get resname]
            puts $fo "$resname$resid $count %"
        }
    }

    close $fo
}



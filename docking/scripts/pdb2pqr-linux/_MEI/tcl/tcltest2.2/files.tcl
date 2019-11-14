# files.tcl -- manage test suite files.

# $Id$

namespace eval tcltest::files {
    array set failFiles {}
}

# tcltest::files::add --
#
#	Add a file to the list of files that have been processed by
#	the test suite.
#
# Arguments:
#	filename - name of file to add to list.
#
# Side Effects:
#	File name list includes filename.
#
# Results:
#	None.

proc tcltest::files::add {filename} {
    variable filelist
    lappend filelist $filename
}

# tcltest::files::getlist --
#
#	Get list of files the test suite has been through.
#
# Arguments:
#	None.
#
# Side Effects:
#	None.
#
# Results:
#	List of files.

proc tcltest::files::getlist {} {
    variable filelist
    return $filelist
}


# tcltest::files::failed --
#
#	Indicate that a particular file has failed tests.
#
# Arguments:
#	filename.
#
# Side Effects:
#	None.
#
# Results:
#	None.

proc tcltest::files::failed {filename} {
    variable failFiles
    set failFiles($filename) 1
}


# tcltest::files::setcreated --
#
#	This is used to indicate that the particular test file has
#	created some files.
#
# Arguments:
#	None.
#
# Side Effects:
#	None.
#
# Results:
#	None.

proc tcltest::files::setcreated { filename {files ""} } {
    variable createdNewFiles
    if { [llength $files] == 0 } {
	return $createdNewFiles($filename)
    }
    set createdNewFiles($filename) $files
}


# tcltest::files::allcreated --
#
#	Returns a list of all the files that have been created.
#
# Arguments:
#	None.
#
# Side Effects:
#	None.
#
# Results:
#	List of files.

proc tcltest::files::allcreated {} {
    variable createdNewFiles
    return [array names createdNewFiles]
}


# tcltest::files::clearcreated --
#
#	Clear list of created files.
#
# Arguments:
#	None.
#
# Side Effects:
#	List of created files is set to empty list.
#
# Results:
#	None.

proc tcltest::files::clearcreated {} {
    variable createdNewFiles
    array unset createdNewFiles
    array set createdNewFiles {}
}


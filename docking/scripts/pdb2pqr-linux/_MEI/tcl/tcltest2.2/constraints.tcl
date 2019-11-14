# constraints.tcl --

# Interface for constraints.

# $Id$

# tcltest::constraints::exists --
#
#	Check to see whether a given constraint exists.
#
# Arguments:
#	constraint.
#
# Side Effects:
#	None.
#
# Results:
#	1 if constraint exists, 0 if it does not.

proc tcltest::constraints::exists {constraint} {
    return [info exists vars::$constraint]
}


# tcltest::constraints::cset --
#
#	Set constraint or check its value.
#
# Arguments:
#	constraint - constraint to set or check.
#	value - optional argument.
#
# Side Effects:
#	Sets constraint if value is given.
#
# Results:
#	None.

proc tcltest::constraints::cset {args} {
    set constraint [lindex $args 0]
    if { [llength $args] == 1 } {
	if { ! [info exists vars::$constraint] } {
	    return 0
	} else {
	    return [set vars::$constraint]
	}
    } else {
	set vars::$constraint [lindex $args 1]
    }
}


proc tcltest::constraints::initconst {constraint} {
    set retval 0
    if { [catch {
	set retval [tcltest::testConstraint $constraint \
			[eval [ConstraintInitializer $constraint]]]
    } err] } {
	puts "DIO CAGNOLINO $err"
    }

    return $retval
}


# tcltest::constraints::getlist --
#
#	Gets a list of all constraints.
#
# Arguments:
#	None.
#
# Side Effects:
#	None.
#
# Results:
#	List of all constraints.

proc tcltest::constraints::getlist {} {
    set reslist {}
    foreach v [info vars vars::*] {
	lappend reslist [namespace tail $v]
    }
    return $reslist
}


# tcltest::constraints::incrskippedbecause --
#
#	Increments the variable used to track how many tests were
#       skipped because of a particular constraint.
#
# Arguments:
#	constraint     The name of the constraint to be modified
#
# Results:
#	Modifies tcltest::skippedBecause; sets the variable to 1 if
#       didn't previously exist - otherwise, it just increments it.
#
# Side effects:
#	None.

proc tcltest::constraints::incrskippedbecause { constraint {value 1} } {
    variable skippedBecause

    if {[info exists skippedBecause($constraint)]} {
	incr skippedBecause($constraint) $value
    } else {
	set skippedBecause($constraint) $value
    }
    return
}


# tcltest::constraints::skippedlist --
#
#	Get list of all constraints that kept tests from running..
#
# Arguments:
#	None.
#
# Side Effects:
#	None.
#
# Results:
#	A list of constraints.

proc tcltest::constraints::skippedlist {} {
    variable skippedBecause
    return [array names skippedBecause]
}


# tcltest::constraints::getskipped --
#
#	Gets number of tests skipped because of a particular
#	constraint.
#
# Arguments:
#	constraint - constraint.
#
# Side Effects:
#	None.
#
# Results:
#	Integer number of tests skipped.

proc tcltest::constraints::getskipped { constraint } {
    variable skippedBecause
    return $skippedBecause($constraint)
}


# tcltest::constraints::clearskippedlist --
#
#	Clears the list of skipped constraints.
#
# Arguments:
#	None.
#
# Side Effects:
#	Resets the list of skipped constraints.
#
# Results:
#	None.

proc tcltest::constraints::clearskippedlist {} {
    variable skippedBecause
    array unset skippedBecause
    array set skippedBecause {}
}


# tcltest::constraints::checktest --
#
#	Check test to see if the constraints are satisfied.  Note that
#	'constraintsvar' has to use upvar to reference the real
#	variable, because these checks actually change the
#	constraints.  Something to fix in the future if possible.
#
# Arguments:
#	name - test name.
#	constraintsvar - constraint to check against.
#
# Side Effects:
#	None.
#
# Results:
#	None.

proc tcltest::constraints::checktest {name constraintsvar} {
    upvar $constraintsvar constraints
    set doTest 0

    # I don't agree with this.  I think that a constraint should
    # either be an artificial construct such as unix || pc, OR it
    # should be a plain old Tcl expression, possibly to be evaluated
    # in its own namespace.  FIXME at some later date when we can toss
    # this stuff out. -davidw

    if {[string match {*[$\[]*} $constraints] != 0} {
	# full expression, e.g. {$foo > [info tclversion]}
	catch {set doTest [uplevel \#0 expr $constraints]}
    } elseif {[regexp {[^.a-zA-Z0-9 \n\r\t]+} $constraints] != 0} {
	# something like {a || b} should be turned into
	# $testConstraints(a) || $testConstraints(b).

	regsub -all {[.\w]+} $constraints {$&} c
	catch {set doTest [namespace eval vars [list expr $c]]}
    } elseif {![catch {llength $constraints}]} {
	# just simple constraints such as {unixOnly fonts}.
	set doTest 1
	foreach constraint $constraints {
	    if { ! [cset $constraint] } {
		set doTest 0
		# store the constraint that kept the test from
		# running
		set constraints $constraint
		break
	    }
	}
    }

    # Return the opposite of doTest
    return [expr {$doTest ? 0 : 1}]
}


# tcltest::constraints::ConstraintInitializer --
#
#	Get or set a script that when evaluated in the tcltest namespace
#	will return a boolean value with which to initialize the
#	associated constraint.
#
# Arguments:
#	constraint - name of the constraint initialized by the script
#	script - the initializer script
#
# Results
#	boolean value of the constraint - enabled or disabled
#
# Side effects:
#	Constraint is initialized for future reference by [test]

proc tcltest::constraints::ConstraintInitializer {constraint {script ""}} {
    variable ConstraintInitializer

    # Check for boolean values
    if {![info complete $script]} {
	return -code error "ConstraintInitializer must be complete script"
    }
    set retval [namespace eval ::tcltest $script]
    cset $constraint $retval

}

# tcltest::constraints::DefineConstraintInitializers --
#
#	Set up the initial constraints (such as unix, pc, and so on).
#
# Arguments:
#	None.
#
# Side Effects:
#	Creates a number of constraints.
#
# Results:
#	None.

proc tcltest::constraints::DefineConstraintInitializers {} {
    ConstraintInitializer singleTestInterp {tcltest::singleProcess}

    # All the 'pc' constraints are here for backward compatibility and
    # are not documented.  They have been replaced with equivalent 'win'
    # constraints.

    ConstraintInitializer unixOnly \
	    {string equal $::tcl_platform(platform) unix}
    ConstraintInitializer macOnly \
	    {string equal $::tcl_platform(platform) macintosh}
    ConstraintInitializer pcOnly \
	    {string equal $::tcl_platform(platform) windows}
    ConstraintInitializer winOnly \
	    {string equal $::tcl_platform(platform) windows}

    ConstraintInitializer unix {tcltest::testConstraint unixOnly}
    ConstraintInitializer mac {tcltest::testConstraint macOnly}
    ConstraintInitializer pc {tcltest::testConstraint pcOnly}
    ConstraintInitializer win {tcltest::testConstraint winOnly}

    ConstraintInitializer unixOrPc \
	    {expr {[tcltest::testConstraint unix] || [tcltest::testConstraint pc]}}
    ConstraintInitializer macOrPc \
	    {expr {[tcltest::testConstraint mac] || [tcltest::testConstraint pc]}}
    ConstraintInitializer unixOrWin \
	    {expr {[tcltest::testConstraint unix] || [tcltest::testConstraint win]}}
    ConstraintInitializer macOrWin \
	    {expr {[tcltest::testConstraint mac] || [tcltest::testConstraint win]}}
    ConstraintInitializer macOrUnix \
	    {expr {[tcltest::testConstraint mac] || [tcltest::testConstraint unix]}}

    ConstraintInitializer nt {string equal $::tcl_platform(os) "Windows NT"}
    ConstraintInitializer 95 {string equal $::tcl_platform(os) "Windows 95"}
    ConstraintInitializer 98 {string equal $::tcl_platform(os) "Windows 98"}

    # The following Constraints switches are used to mark tests that
    # should work, but have been temporarily disabled on certain
    # platforms because they don't and we haven't gotten around to
    # fixing the underlying problem.

    ConstraintInitializer tempNotPc {expr {![tcltest::testConstraint pc]}}
    ConstraintInitializer tempNotWin {expr {![tcltest::testConstraint win]}}
    ConstraintInitializer tempNotMac {expr {![tcltest::testConstraint mac]}}
    ConstraintInitializer tempNotUnix {expr {![tcltest::testConstraint unix]}}

    # The following Constraints switches are used to mark tests that
    # crash on certain platforms, so that they can be reactivated again
    # when the underlying problem is fixed.

    ConstraintInitializer pcCrash {expr {![tcltest::testConstraint pc]}}
    ConstraintInitializer winCrash {expr {![tcltest::testConstraint win]}}
    ConstraintInitializer macCrash {expr {![tcltest::testConstraint mac]}}
    ConstraintInitializer unixCrash {expr {![tcltest::testConstraint unix]}}

    # Skip empty tests

    ConstraintInitializer emptyTest {format 0}

    # By default, tests that expose known bugs are skipped.

    ConstraintInitializer knownBug {format 0}

    # By default, non-portable tests are skipped.

    ConstraintInitializer nonPortable {format 0}

    # Some tests require user interaction.

    ConstraintInitializer userInteraction {format 0}

    # Some tests must be skipped if the interpreter is not in
    # interactive mode

    ConstraintInitializer interactive \
	    {expr {[info exists ::tcl_interactive] && $::tcl_interactive}}

    # Some tests can only be run if the installation came from a CD
    # image instead of a web image.  Some tests must be skipped if you
    # are running as root on Unix.  Other tests can only be run if you
    # are running as root on Unix.

    ConstraintInitializer root {expr \
	    {[string equal unix $::tcl_platform(platform)]
	    && ([string equal root $::tcl_platform(user)]
		|| [string equal "" $::tcl_platform(user)])}}
    ConstraintInitializer notRoot {expr {![tcltest::testConstraint root]}}

    # Set nonBlockFiles constraint: 1 means this platform supports
    # setting files into nonblocking mode.

    ConstraintInitializer nonBlockFiles {
	    set code [expr {[catch {set f [open defs r]}] 
		    || [catch {fconfigure $f -blocking off}]}]
	    catch {close $f}
	    set code
    }

    # Set asyncPipeClose constraint: 1 means this platform supports
    # async flush and async close on a pipe.
    #
    # Test for SCO Unix - cannot run async flushing tests because a
    # potential problem with select is apparently interfering.
    # (Mark Diekhans).

    ConstraintInitializer asyncPipeClose {expr {
	    !([string equal unix $::tcl_platform(platform)] 
	    && ([catch {exec uname -X | fgrep {Release = 3.2v}}] == 0))}}

    # Test to see if we have a broken version of sprintf with respect
    # to the "e" format of floating-point numbers.

    ConstraintInitializer eformat {string equal [format %g 5e-5] 5e-05}

    # Test to see if execed commands such as cat, echo, rm and so forth
    # are present on this machine.

    ConstraintInitializer unixExecs {
	set code 1
        if {[string equal macintosh $::tcl_platform(platform)]} {
	    set code 0
        }
        if {[string equal windows $::tcl_platform(platform)]} {
	    if {[catch {
	        set file _tcl_test_remove_me.txt
	        makeFile {hello} $file
	    }]} {
	        set code 0
	    } elseif {
	        [catch {exec cat $file}] ||
	        [catch {exec echo hello}] ||
	        [catch {exec sh -c echo hello}] ||
	        [catch {exec wc $file}] ||
	        [catch {exec sleep 1}] ||
	        [catch {exec echo abc > $file}] ||
	        [catch {exec chmod 644 $file}] ||
	        [catch {exec rm $file}] ||
	        [llength [auto_execok mkdir]] == 0 ||
	        [llength [auto_execok fgrep]] == 0 ||
	        [llength [auto_execok grep]] == 0 ||
	        [llength [auto_execok ps]] == 0
	    } {
	        set code 0
	    }
	    removeFile $file
        }
	set code
    }

    ConstraintInitializer stdio {
	set code 0
	if {![catch {set f [open "|[list [interpreter]]" w]}]} {
	    if {![catch {puts $f exit}]} {
		if {![catch {close $f}]} {
		    set code 1
		}
	    }
	}
	set code
    }

    # Deliberately call socket with the wrong number of arguments.  The
    # error message you get will indicate whether sockets are available
    # on this system.

    ConstraintInitializer socket {
	catch {socket} msg
	string compare $msg "sockets are not available on this system"
    }

    # Check for internationalization
    ConstraintInitializer hasIsoLocale {
	if {[llength [info commands testlocale]] == 0} {
	    set code 0
	} else {
	    set code [string length [SetIso8859_1_Locale]]
	    RestoreLocale
	}
	set code
    }

}

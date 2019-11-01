namespace eval tcltest::testresults {

    # The 'TestNames' variable doesn't work for us because it's an
    # array.
    set testNames [list]

    array set ResID {}

    # Since this is a pretty long list, we keep it in one place
    # here, and loop through it when needed.
    set testAPIvars {
	testDescription
	testPassFail
	testBody
	testDontRun

	testSkipped

	testSetup
	testScriptFailure
	testScriptMatch
	testActualAnswer
	testMatch
	testResult
	testMsg
	testReturnCode
	testExpectedReturnCodes
	testErrorInfo
	testErrorCode

	testCodeFailure
	testOutputFailure
	testErrorFailure

	testOutputMatch
	testOutData
	testOutput

	testErrorMatch
	testErrorData
	testErrorOutput

	testCleanupMsg
	testCore
    }
    foreach var $testAPIvars {
	array set $var {}
    }

    # initialize numTests array to keep track of the number of tests
    # that pass, fail, and are skipped.
    array set numTests [list Total 0 Passed 0 Skipped 0 Failed 0]
}

# tcltest::testresults::newresult --
#
#	Start recording information about a new result.  Must be
#	called prior to tcltest::result.
#
# Arguments:
#	name - test name.
#
# Side Effects:
#	Sets the ResID variable, which acts as a unique ID for each
#	result.
#
# Results:
#	None.

proc tcltest::testresults::newresult {name} {
    variable testNames
    variable ResID

    lappend testNames $name
    set ResID($name) [expr {[llength $testNames] - 1}]
}

# tcltest::testresults::result --
#
#	Store information about a particular test.
#
# Arguments:
#	var - variable name to set.
#	name - test name.
#	data - data to store.
#
# Side Effects:
#	None.
#
# Results:
#	None.

proc tcltest::testresults::result {var name data} {
    variable testAPIvars
    variable ResID

    foreach tstvar $testAPIvars {
	variable $tstvar
    }

    set [set var]($ResID($name)) $data
}


# tcltest::testresults::results --
#
#	Takes a subcommand as the first argument, and based on that,
#	returns some information about test results.
#
# Arguments:
#	cmd - one of tests vars, exists, get or clear
#	args - additional arguments, such as the test name, the
#	variable to query, and so on.  See the tcltest man page.
#
# Side Effects:
#	The 'clear' command eliminates the results of the specified test.
#
# Results:
#	tests - returns information on all tests.
#	vars - returns the variables available for a given test.
#	exists - reports whether a variable exists for a given test.
#	get - fetches the information about a variable for a given test.
#	clear - clears all information about a given test, and removes
#	it from the list of tests.

proc tcltest::testresults::results {cmd args} {
    variable testNames
    variable testAPIvars
    variable numTests

    foreach var $testAPIvars {
	variable $var
    }

    if { $cmd == "tests" } {
	return $testNames
    }

    set res {}
    set testname [lindex $args 0]
    # Fetch the list index of the test name.
    set id [lsearch -all $testNames $testname]
    # If there was more than one result, either one is specified or we
    # error out.

    if { [llength $id] > 1 } {
	set num [lindex $args 1]
	if { ![string is int $num] || $num == "" } {
	    error "Multiple tests correspond - please specify one: 0-[expr [llength $id] - 1]"
	}
	set id [lindex $id $num]
    } elseif {$id < 0} {
	error "No $testname test recorded"
    }

    switch -exact $cmd {
	vars {
	    foreach var $testAPIvars {
		if { [info exists [set var]($id)] } {
		    lappend res $var
		}
	    }
	    return $res
	}

	exists {
	    set varname [lindex $args end]
	    return [info exists [set varname]($id)]
	}

	get {
	    set varname [lindex $args end]
	    return [set [set varname]($id)]
	}

	clear {
	    foreach var $testAPIvars {
		if { [info exists [set var]($id)] } {
		    unset [set var]($id)
		}
	    }
	    set testNames [lreplace $testNames $id $id]
	    return
	}

	default {
	    error "bad option \"$cmd\": must be tests, vars, exists, get, or clear"
	}
    }
}


# tcltest::testresults::tally --
#
#	Return information about the number of tests that have passed,
#	failed, been skipped, and the total.
#
# Arguments:
#	type - type of tally to return: passed, failed, skipped, total.
#
# Side Effects:
#	None.
#
# Results:
#	A number corresponding to the type of information sought.

proc tcltest::testresults::tally {type} {
    variable numTests
    return $numTests([string totitle $type])
}


# tcltest::testresults::incrtotal --
#
#	Increase the count for the number of tests that have either
#	Passed, Failed, Skipped, or the Total number of tests.
#
# Arguments:
#	None.
#
# Side Effects:
#	None.
#
# Results:
#	None.

proc tcltest::testresults::incrtotal {type {num 1}} {
    variable numTests
    incr numTests($type) $num
}
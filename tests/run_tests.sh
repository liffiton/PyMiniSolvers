#!/bin/bash
#
# run_tests.sh -- Run regression tests
#
# Author: Mark Liffiton
# Date: June 2011
#

mode="run"
TESTSFILE="all_tests.txt"

# setup indicator characters
chrPass="[32m*[0m"
chrSortSame="[33m^[0m"
chrStdErr="[34mo[0m"

if [ "$1" != "" ] ; then
    mode=$1
fi

if [ "$2" != "" ] && [ -f $2 ] ; then
    TESTSFILE="$2"
fi

checkfiles() {
    f1=$1
    f2=$2

    # first compare sizes
    size1=`du -b $f1 | cut -f 1`
    size2=`du -b $f2 | cut -f 1`
    if [ "$size1" -ne "$size2" ] ; then
        echo
        echo "  [31mOutputs differ (size).[0m"
        return 11
    fi

    # then check outputs directly, and if still different then sort and check again
    cmp -s $f1 $f2
    if [ $? -ne 0 ] ; then
        cmp -s <(sort $f1) <(sort $f2)
        if [ $? -ne 0 ] ; then
            echo
            echo "  [31mOutputs differ (contents).[0m"
            return 12
        else
            # outputs not equivalent, but sort to same contents
            return 1
        fi
    fi

    # everything checks out
    return 0
}

viewdiff() {
    f1=$1
    f2=$2
    read -p "  View diff? (T for terminal, V for vimdiff) " -n 1
    echo
    if [[ $REPLY =~ ^[Tt]$ ]]; then
        $DIFF -C3 $f1 $f2
    elif [[ $REPLY =~ ^[Vv]$ ]]; then
        vimdiff $f1 $f2
    fi
}

updateout() {
    f1=$1
    f2=$2
    read -p "  Store new output as correct? " -n 1
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "  [33mmv $f2 $f1[0m"
        mv $f2 $f1
    fi
}

runtest() {
    cmd=$1
    outfile=$2

    [ $verbose -eq 1 ] && echo -e "\nRunning test: $cmd > $outfile.NEW" || echo -n "."

    [ $verbose -eq 1 ] && $cmd > $outfile.NEW || $cmd > $outfile.NEW 2> $outfile.ERR

    if [ -s $outfile.ERR ] ; then
        [ $verbose -eq 0 ] && echo -ne "\b$chrStdErr."  # note that stderr produced output
        (( stdErrOutput++ ))
    fi
    rm -f $outfile.ERR

    if [ $? -gt 128 ] ; then
        # Fatal error! (likely a segfault)
        echo "  [37;41mTest failed (error code returned):[0m $cmd"
    else
        checkfiles $outfile $outfile.NEW
        ret=$?

        if [ $ret -eq 0 ] ; then
            [ $verbose -eq 1 ] && echo "  [32mTest passed.[0m" || echo -ne "\b$chrPass"  # \b = backspace, to overwrite the .
            (( passed++ ))
        elif [ $ret -eq 1 ] ; then
            [ $verbose -eq 1 ] && echo -e "\n  [33mOutputs not equivalent, but sort to same contents.[0m" || echo -ne "\b$chrSortSame"
            (( passed++ ))
            (( sortedSame++ ))
        else
            echo
            echo "  [37;41mTest failed:[0m $cmd"
            viewdiff $outfile $outfile.NEW
            updateout $outfile $outfile.NEW
        fi
    fi

    rm -f $outfile.NEW

    (( total++ ))
}

generate() {
    cmd=$1
    outfile=$2
    echo "Generating: $cmd > [32m$outfile[0m"
    $cmd > $outfile
}

##################################################################################

# Make sure we're in the tests directory
cd `dirname $0`

# Test for existence of our test file
if [ ! -e "$TESTSFILE" ] ; then
    echo "ERROR: Tests file $TESTSFILE does not exist."
    exit 1
fi

# Check for colordiff
[ -x "`which colordiff 2>/dev/null`" ] && DIFF="colordiff" || DIFF="diff"

# Set our mode (run as default) and verbosity (0 as default)
verbose=0
case "$mode" in
    run | "")
        mode="run"
        echo "Running all tests from $TESTSFILE."
        ;;
    runverbose)
        mode="run"
        verbose=1
        echo "Running all tests from $TESTSFILE."
        ;;
    regenerate)
        read -p "Are you sure you want to regenerate all test outputs (y/n)? "
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Exiting."
            exit 1
        fi
        mode="regenerate" ;;
    nocheck)
        mode="nocheck"
        echo "Running all tests from $TESTSFILE (skipping results checks)"
        ;;
    *)
        echo "Invalid mode: $mode"
        echo "Options: run, runverbose, regenerate, nocheck"
        exit 1
        ;;
esac

total=0
passed=0
stdErrOutput=0
sortedSame=0

# Run a test for each line in our tests file
# Use a separate file descriptor so we still have STDIN for asking questions via read
exec 4< $TESTSFILE  # FD4 attaches to $TESTSFILE
while read LINE <&4; do
    # Skip blank lines
    if [ "$LINE" = "" ] ; then
        continue
    fi
    # Skip "commented" lines
    if [ "${LINE:0:1}" = "#" ] ; then
        continue
    fi

    # Split the line into an array
    IFS=, read cmd outfile <<< "$LINE"

    exe=`cut -d " " -f 1 <<< "$cmd"`
    if [ ! -x $exe ] ; then
        echo "ERROR: $exe is not an executable file.  Do you need to run make?"
        exit 1
    fi

    if [ "$mode" = "run" ] ; then
        runtest "$cmd" "$outfile"
    elif [ "$mode" = "regenerate" ] ; then
        generate "$cmd" "$outfile"
    elif [ "$mode" = "nocheck" ] ; then
        $cmd > /dev/null
        echo "."
    fi

done

if [ "$mode" = "run" ] ; then
    echo
    echo " $chrPass : $passed / $total   Passed"
    if [ $sortedSame -gt 0 ] ; then
        echo " $chrSortSame : $sortedSame   [Different order, but same output]"
    fi
    if [ $stdErrOutput -gt 0 ] ; then
        echo " $chrStdErr : $stdErrOutput   [Produced output to STDERR]"
    fi
    echo
fi

exit 0


#!/bin/sh

errorout()
{
    echo $1 #> /dev/stderr
    exit 1
}

[ $# -ne 3 ] || errorout "Please supply exactly two arguments to the test script [BINARY DIR] [TEMP DIR]."

BINDIR=$1
TMPDIR=$2

${BINDIR}/sviper -h
[ $? -eq 0 ] || errorout "Could not print the help page"

exit 0

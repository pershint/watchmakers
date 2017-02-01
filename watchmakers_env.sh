#!/bin/bash

#echo $BASH_SOURCE
LOC=`dirname $BASH_SOURCE`
#LOCATION="$PWD/$LOC/watchmakers.py "
LOCATION="$LOC/watchmakers.py "
#echo $LOCATION
alias watch="python $LOCATION"

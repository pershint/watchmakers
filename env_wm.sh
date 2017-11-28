#!/bin/sh
WATCHENV=/home/liz/watchmakers
PATH=$WATCHENV:$PATH
PYTHONPATH=$WATCHENV:$PYTHONPATH
export WATCHENV PATH PYTHONPATH
alias watch="python $WATCHENV/watchmakers.py "

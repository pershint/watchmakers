#!/bin/csh
setenv WATCHENV /home/liz/watchmakers
setenv PATH "$WATCHENV:$PATH"

if ({$?PYTHONPATH}) then
  setenv PYTHONPATH "$WATCHENV:$PYTHONPATH"
else
  setenv PYTHONPATH "$WATCHENV"
endif
alias watch "python /watchmakers.py "

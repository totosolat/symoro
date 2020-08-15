#!/bin/sh

DIR="$(dirname $0)"

PYTHONPATH="$DIR":"$PYTHONPATH" python2 "$DIR"/bin/symoro-bin.py


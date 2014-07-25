#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
test_out=$DIR/test_out
rm -rf $test_out
$DIR/../bin/polymarker.rb -c $DIR/data/PST130_7067.fa -s  $DIR/data/PST130_7067.csv -r $DIR/data/PST130_7067.fa  -o $test_out -e affine:local -a arm_selection_morex

#!/bin/bash

time sort -S 100G -r -g -k5 --parallel=10 -T /localscratch/tmp $1 > $1.sorted

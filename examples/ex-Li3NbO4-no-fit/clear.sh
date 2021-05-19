#!/bin/bash

find . -name '*out' -exec rm -rv {} \;
find . -name 'verlet.*' -exec rm -v {} \;

#!/bin/bash

HERE=$(pwd)
OSPREY_bin=/app/misc/OSPREY/bin
java -cp $OSPREY_bin -Xmx4G KStar -t 3 -c $HERE/KStar.cfg doDEE $1 $2

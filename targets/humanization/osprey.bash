#!/bin/bash

target="cetuximab"
java -Xmx4000M KStar -t 3 -c KStar.cfg doDEE $target.System.cfg $target.DEE.cfg > $target.log

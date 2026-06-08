#!/bin/bash

target="lst_cwbd"
java -Xmx4000M KStar -t 3 -c KStar.cfg doDEE $target.System.cfg $target.DEE.cfg > $target"_DEE.log"

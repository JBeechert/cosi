#!/bin/bash
for number in {1..50}
do
mimrec -c mimrec_extract_onGP_cuts_10dets.cfg -f ../AllSky.Al26.160517-160605.p1.inc${number}.id1.tra.gz -x -n -o AllSky.Al26.160517-160605.p1.inc${number}.id1_onGPcuts.tra.gz
sleep 3
done
exit 0

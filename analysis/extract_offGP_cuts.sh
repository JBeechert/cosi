#!/bin/bash
for number in {1..50}
do
mimrec -c mimrec_extract_offGP_cuts_lat0_90_long-180_80_10dets.cfg -f ../AllSky.Al26.160517-160605.p1.inc${number}.id1.tra.gz -x -n -o AllSky.Al26.160517-160605.p1.inc${number}.id1_offGPcuts_lat0_90_long-180_80.tra.gz
mimrec -c mimrec_extract_offGP_cuts_lat-85_5_long0_30_10dets.cfg -f ../AllSky.Al26.160517-160605.p1.inc${number}.id1.tra.gz -x -n -o AllSky.Al26.160517-160605.p1.inc${number}.id1_offGPcuts_lat-85_5_long0_30.tra.gz
mimrec -c mimrec_extract_offGP_cuts_lat85_5_long0_30_10dets.cfg -f ../AllSky.Al26.160517-160605.p1.inc${number}.id1.tra.gz -x -n -o AllSky.Al26.160517-160605.p1.inc${number}.id1_offGPcuts_lat85_5_long0_30.tra.gz

done
exit 0

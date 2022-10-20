#!/bin/bash

labelsize="100"

hopsize="3"

diameter="3"

queryNum="1"

pathlength="3"

querysize="8"

#mode 0:fast test 1:homo 2:strong simulation
mode="2"


##### Original Test Script without BloomFilter+SGX

#rm -rf results/*

rm CMakeCache.txt

cmake .

make

for((i=1;i<=$queryNum;i++));
do
./main query/Query-Slashdot/S$i-$labelsize graph/slashdot$labelsize $labelsize $hopsize $pathlength $mode
done

for((i=11;i<=$queryNum;i++));
do
./main query/Query-DBLP/D$i-$labelsize graph/dblp-un$labelsize $labelsize $hopsize $pathlength $mode
done

for((i=1;i<=$queryNum;i++));
do
./main query/Query-Twitter/T$i-$labelsize graph/twitter$labelsize $labelsize $hopsize $pathlength $mode
done







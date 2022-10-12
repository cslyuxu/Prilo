labelsize="100"

hopsize="3"

diameter="4"

queryNum="10"

pathlength="3"

querysize="8"

threshold="5"

## To generate queries for Slashdot, DBLP and Twitter, use PPMatch/PPMatch/QGen.sh

#mode 0:fast test 1:homo 2:strong simulation
mode="0"

make clean

make

for((i=1;i<=$queryNum;i++));
do
./PPMatchApp PPMatch/query/Query-Slashdot/S$i-$labelsize PPMatch/graph/slashdot$labelsize $labelsize $hopsize $pathlength $mode $threshold
done

for((i=11;i<=$queryNum;i++));
do
./PPMatchApp PPMatch/query/Query-DBLP/D$i-$labelsize PPMatch/graph/dblp-un$labelsize $labelsize $hopsize $pathlength $mode $threshold
done

for((i=1;i<=$queryNum;i++));
do
./PPMatchApp PPMatch/query/Query-Twitter/T$i-$labelsize PPMatch/graph/twitter$labelsize $labelsize $hopsize $pathlength $mode $threshold
done

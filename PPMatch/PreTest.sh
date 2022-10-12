#!/bin/bash

labelsize="100"

hopsize="3"

diameter="3"

queryNum="1"

pathlength="3"

querysize="8"

#mode 0:fast test 1:homo 2:strong simulation
mode="2"


##### GrapghTransformation



#cd GraphTransformation

#./test.sh $labelsize

#cd ../

#cp -f GraphTransformation/graph/slashdot$labelsize GraphIndex/graph/slashdot$labelsize
#cp -f GraphTransformation/graph/dblp-un$labelsize GraphIndex/graph/dblp-un$labelsize
#cp -f GraphTransformation/graph/twitter$labelsize GraphIndex/graph/twitter$labelsize
#cp -f GraphTransformation/graph/slashdot$labelsize QueryGenerator/graph/slashdot$labelsize
#cp -f GraphTransformation/graph/dblp-un$labelsize QueryGenerator/graph/dblp-un$labelsize
#cp -f GraphTransformation/graph/twitter$labelsize QueryGenerator/graph/twitter$labelsize
#cp -f GraphTransformation/graph/slashdot$labelsize graph/slashdot$labelsize
#cp -f GraphTransformation/graph/dblp-un$labelsize graph/dblp-un$labelsize
#cp -f GraphTransformation/graph/twitter$labelsize graph/twitter$labelsize



##### QueryGenerator

#cd QueryGenerator

#rm -rf Query-Slashdot/*

#rm -rf Query-DBLP/*

#rm -rf Query-Twitter/*

#./test.sh $diameter $queryNum $querysize $labelsize

#cd ../

#rm -rf query/Query-Slashdot/*

#rm -rf query/Query-DBLP/*

#rm -rf query/Query-Twitter/*


#cp -rf QueryGenerator/Query-Slashdot/* PathGenerator/Query-Slashdot/

#cp -rf QueryGenerator/Query-DBLP/* PathGenerator/Query-DBLP/

#cp -rf QueryGenerator/Query-Twitter/* PathGenerator/Query-Twitter/

#cp -rf QueryGenerator/Query-Slashdot/* query/Query-Slashdot/

#cp -rf QueryGenerator/Query-DBLP/* query/Query-DBLP/

#cp -rf QueryGenerator/Query-Twitter/* query/Query-Twitter/



##### Test without BloomFilter+SGX

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







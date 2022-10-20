#!/bin/bash

labelsize="64"

diameter="4"

queryNum="10"

querysize="8"


##### GrapghTransformation

cd GraphTransformation

./test.sh $labelsize

cd ../

cp -f GraphTransformation/graph/slashdot$labelsize QueryGenerator/graph/slashdot$labelsize
cp -f GraphTransformation/graph/dblp-un$labelsize QueryGenerator/graph/dblp-un$labelsize
cp -f GraphTransformation/graph/twitter$labelsize QueryGenerator/graph/twitter$labelsize
cp -f GraphTransformation/graph/slashdot$labelsize graph/slashdot$labelsize
cp -f GraphTransformation/graph/dblp-un$labelsize graph/dblp-un$labelsize
cp -f GraphTransformation/graph/twitter$labelsize graph/twitter$labelsize



##### QueryGenerator

cd QueryGenerator

rm -rf Query-Slashdot/*

rm -rf Query-DBLP/*

rm -rf Query-Twitter/*

./test.sh $diameter $queryNum $querysize $labelsize

cd ../

rm -rf query/Query-Slashdot/*

rm -rf query/Query-DBLP/*

rm -rf query/Query-Twitter/*

cp -rf QueryGenerator/Query-Slashdot/* query/Query-Slashdot/

cp -rf QueryGenerator/Query-DBLP/* query/Query-DBLP/

cp -rf QueryGenerator/Query-Twitter/* query/Query-Twitter/












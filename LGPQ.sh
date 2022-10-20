#Label size          Ssim: 64      Hom: Slashdot---100, DBLP---150, Twitter---100
labelsize="64" 

queryNum="10"

#h for Twiglet: 3-5 
twigletlength="3"

#t for TreeBF: 5, 15, 25
threshold="15"

#mode 0: Test PPCR (TreeBF+Twiglet) on Ssim due to more evaluated balls     1: Hom      2: Ssim
mode="0"

make clean

make

for((i=11;i<=$queryNum;i++));
do
./LGPQApp LGPQ/query/Query-Slashdot/S$i-$labelsize LGPQ/graph/slashdot$labelsize $labelsize $twigletlength $mode $threshold
done

for((i=1;i<=$queryNum;i++));
do
./LGPQApp LGPQ/query/Query-DBLP/D$i-$labelsize LGPQ/graph/dblp-un$labelsize $labelsize $twigletlength $mode $threshold
done

for((i=11;i<=$queryNum;i++));
do
./LGPQApp LGPQ/query/Query-Twitter/T$i-$labelsize LGPQ/graph/twitter$labelsize $labelsize $twigletlength $mode $threshold
done

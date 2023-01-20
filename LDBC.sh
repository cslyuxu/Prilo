#Label size 
labelsize="213" 

queryNum="10"

#h for Twiglet: 3-5 
twigletlength="3"

#t for TreeBF: 5, 15, 25
threshold="15"

#mode 0: Test PPCR (TreeBF+Twiglet) on Ssim due to more evaluated balls     1: Hom      2: Ssim
mode="1"

make clean

make

for((i=1;i<=$queryNum;i++));
do
./LGPQApp LGPQ/query/Query-LDBC/L$i LGPQ/graph/LDBC_sf1 $labelsize $twigletlength $mode $threshold
done


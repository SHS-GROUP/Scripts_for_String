#!/bin/sh

IPATH=/home/pengfeil/Projs/SLO_QMMM/6b-string_l546a_l754a/
cycle=4
imgs=19
fc=100.00

for ((a=1;a<=${cycle};a++))
do
b=$(($a-1));

mkdir -p $a
cd $a

# About the file containing force constants
if [ -f force.$a.dat ]
then
rm force.$a.dat
fi

for ((i=1;i<=${imgs};i++))
do
echo "${fc} ${fc} ${fc}" >> force.$a.dat
done

# About the files containing equlibrium distances
cp $IPATH/optstring_calcs$b/newconstr.dat newconstr.$a.dat

# About the files containing actual sampling data
for ((img=1;img<=${imgs};img++))
do
cp $IPATH/optstring_calcs$b/s${b}_i${img}_collect.dat distperframe${img}.dat
done

cd ..
done


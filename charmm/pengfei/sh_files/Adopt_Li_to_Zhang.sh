#!/bin/sh

IPATH=/home/pengfeil/Projs/String_test/4-mastro/transfer_to_maestro/running

for ((a=1;a<=5;a++))
do
b=$(($a-1));

mkdir -p $a
cd $a
cp $IPATH/analysis/force.dat force.$a.dat
cp $IPATH/optstring_calcs$b/newconstr.dat newconstr.$a.dat

for ((img=1;img<=18;img++))
do
mv $IPATH/analysis/s${a}_i${img}_collect.dat distperframe${img}.dat
done

cd ..

done


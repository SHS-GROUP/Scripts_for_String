#!/bin/sh

# You NEED dis_atnums.txt file in the running directory ($pathn) to specify the atom numbers in each reaction coordinate
# AND react_link.psf, react_qmmin.crd, product_qmmin.crd files in the interpol directory under the running path ($pathn)

# Set up the running directories in shscluster and xsede
pathn=/home/pengfeil/Projs/String_test/4-mastro/transfer_to_maestro/running_test
xsedepath=pfliuiuc@comet.sdsc.xsede.org:/oasis/projects/nsf/uic317/pfliuiuc/Projs/String_test/4-mastro/transfer_to_maestro/running_test
winds=18 # Number of images
dim=3 # Number of reaction coordinates
forcecons=100.0000 # Force constants of the harmonic strings
forceamber="100.0"

if [ -f iniconstr.dat ]
then
    rm iniconstr.dat
fi

if [ -f iniconstr.RST ]
then
    rm iniconstr.RST
fi

if [ -f force.dat ]
then
    rm force.dat
fi

# Do the interpolation and generate the initial constants file
for ((a=1;a<=${winds};a++))
do
if [ -d f$a ]
then
    rm -r f$a
fi
mkdir -p f$a
cd f$a
cp ../*.crd .
cp ../react_link.psf .
crd_inter.py -p react_link.psf --rc=react_qmmin.crd --pc=product_qmmin.crd -o interpol_$a.crd --prog=charmm --img=18 --wind=$a
cal_dis.py -i $pathn/dis_atnums.txt -p react_link.psf -c interpol_$a.crd -o interpol_dis_$a.dat
awk '{print $3, ORS=" "} END {print "\n"}' interpol_dis_$a.dat >> $pathn/interpol/iniconstr.dat
awk '{printf("&rst iat=%d,%d, r1=0., r2=%7.4f, r3=%7.4f, r4=100., rk2='${forceamber}', rk3='${forceamber}',/\n", $1, $2, $3, $3)}' interpol_dis_$a.dat >> $pathn/interpol/iniconstr.RST
#cp ../*.crd .
#cp ../reactant.psf .
#crd_inter.py -p reactant.psf --rc=reactant.crd --pc=product.crd -o interpol_$a.crd --prog=charmm --img=18 --wind=$a
#cal_dis.py -i $pathn/dis_atnums.txt -p reactant.psf -c interpol_$a.crd -o interpol_dis_$a.dat
#awk '{print $3, ORS=" "} END {print "\n"}' interpol_dis_$a.dat >> $pathn/interpol/iniconstr.dat
cd ..
done

# Generate the force constant file
for ((a=1;a<=${winds};a++))
do
for ((b=1;b<=${dim};b++))
do
echo -n "   $forcecons"  >> force.dat
done
echo " " >> force.dat
done

# Copy the interpol dir back to xsede
scp -r $pathn/interpol $xsedepath/


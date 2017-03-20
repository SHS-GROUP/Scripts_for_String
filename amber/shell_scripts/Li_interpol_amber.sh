#!/bin/sh

# Should run this script in the interpol directory, which is under the running directory ($pathn)
# You NEED dis_atnums.txt file in the running directory ($pathn) to specify the atom numbers in each reaction coordinate
# AND react_link.psf, react_qmmin.crd, product_qmmin.crd files in the interpol directory under the running path ($pathn)

# Set up the running directories in shscluster and xsede
pathn="/home/pengfeil/Projs/SLO_QMMM/TESTS/string_running"
xsedepath="pfliuiuc@comet.sdsc.xsede.org:/oasis/projects/nsf/uic317/pfliuiuc/Projs/SLO/QMMM/test/string_running"
imgs=20 # Number of images
dim=3 # Number of reaction coordinates
forcecons=100.0000 # Force constants of the harmonic strings for the k in (1/2)*k(r-req)**2 (which CHARMM has and WHAM code has)
aforcecons=`echo '0.5*'${forcecons}'' | bc`
# Since AMBER has the potential function form as k(r-req)**2 other than (1/2)*k(r-req)**2,
# so the k values in RST files are half of the forcecons
react_top="pre_react.prmtop"
react_crd="pre_react_min.rst"
prod_crd="pre_prod.rst"
atmpairf="string_atompairs.txt"

# Delete the former iniconstr.dat and force.dat files
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
for ((a=1;a<=${imgs};a++))
do

if [ -d f$a ]
then
    rm -rf f$a
fi

mkdir -p f$a
cd f$a
cp $pathn/interpol/$react_top .
cp $pathn/interpol/$react_crd .
cp $pathn/interpol/$prod_crd .
crd_inter.py -p $react_top --rc=$react_crd --pc=$prod_crd -o interpol_$a.inpcrd --prog=amber --img=$imgs --wind=$a
cal_dis.py -i $pathn/$atmpairf -p $react_top -c interpol_$a.inpcrd -o interpol_dis_$a.dat --RST=ini_${a}.RST -k $aforcecons
cat ini_${a}.RST >> ../iniconstr.RST
awk '{print $3, ORS=" "} END {print "\n"}' interpol_dis_$a.dat >> $pathn/interpol/iniconstr.dat
cd ..
done

# Generate the force constant file
for ((a=1;a<=${imgs};a++))
do
for ((b=1;b<=${dim};b++))
do
echo -n "   $forcecons"  >> force.dat
done
echo " " >> force.dat
done

# Copy the interpol dir back to xsede
scp -r $pathn/interpol $xsedepath/


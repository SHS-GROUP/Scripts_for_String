#!/bin/sh

step=$1 # Which step are you in? s0, s1, or s2, etc.?
winds=18  # Windows is the number of images
dim=3 # Dimension is the number of reaction coordinates

# Running directory path in shscluster
pathn=/home/pengfeil/Projs/String_test/4-mastro/transfer_to_maestro/running_test
# Running directory path in xsede
xsedepath=pfliuiuc@comet.sdsc.xsede.org:/oasis/projects/nsf/uic317/pfliuiuc/Projs/String_test/4-mastro/transfer_to_maestro/running_test

scp -r $xsedepath/tmp/ s$step.xsede

if [ -d $pathn/optstring_calcs$step ]
then
    rm -r $pathn/optstring_calcs$step
fi

mkdir -p $pathn/optstring_calcs$step

# Atom numbers of each reaction coordiante are stored in the $pathn/dis_atnums.txt file
for ((a=1;a<=${winds};a++));do
    echo "****Windows${a}****";
    cal_dis.py -i $pathn/dis_atnums.txt -p $pathn/s$step.xsede/react_link.psf -c $pathn/s$step.xsede/f$a/${step}_production.dcd -o $pathn/optstring_calcs$step/$a.dat -t $pathn/optstring_calcs$step/s${step}_i${a}_collect.dat;
done

# Regenerate the string
cd $pathn/optstring_calcs$step;
rege_string.py -i $pathn/dis_atnums.txt -d $dim --img=$winds --nimg=$winds;

# Copy the newconstr.dat back to Xsede
scp $pathn/optstring_calcs$step/newconstr.dat $pathn/optstring_calcs$step/newconstr.RST $xsedepath/s$step/


#!/bin/sh

step=$1 # Which step are you in? s0, s1, or s2, etc.?
winds=19  # Windows is the number of images
dim=3 # Dimension is the number of reaction coordinates

# Running directory path in shscluster
pathn="/home/pengfeil/Projs/SLO_QMMM/6c-string_running"
# Running directory path in xsede
xsedepath="pfliuiuc@comet.sdsc.xsede.org:/oasis/projects/nsf/uic317/pfliuiuc/Projs/SLO/QMMM/4b-string_running"
forcecons=200.0000 # Force constants of the harmonic strings for the k in (1/2)*k(r-req)**2 (which CHARMM has and WHAM code has)
aforcecons=`echo '0.5*'${forcecons}'' | bc`
# Since AMBER has the potential function form as k(r-req)**2 other than (1/2)*k(r-req)**2,
# so the k values in RST files are half of the forcecons
prmtop="slo_3pzw_fe3_lid_react"

if [ "$2" == "total" ]
then
    if [ -d $pathn/s$step.xsede ]
    then
        rm -r $pathn/s$step.xsede
    fi

    scp -r $xsedepath/tmp/ $pathn/s$step.xsede

    if [ -d $pathn/optstring_calcs$step ]
    then
        rm -r $pathn/optstring_calcs$step
    fi

    mkdir -p $pathn/optstring_calcs$step

    #Atom numbers of each reaction coordiante are stored in the $pathn/dis_atnums.txt file
    for ((a=1;a<=${winds};a++));do
        echo "****Windows${a}****";
        cal_dis.py -i $pathn/string_atompairs.txt -p $pathn/s$step.xsede/${prmtop}.prmtop -c $pathn/s$step.xsede/f$a/${step}_production.netcdf -o $pathn/optstring_calcs$step/$a.dat -t $pathn/optstring_calcs$step/s${step}_i${a}_collect.dat;
    done

    #Global binning
    for ((a=1;a<=${dim};a++))
    do
        glo_bin.py -c $1 -w ${winds} -b 20 -d ${a}
    done

    # Regenerate the string
    cd $pathn/optstring_calcs$step;
    rege_string.py -i $pathn/string_atompairs.txt -d ${dim} --img=${winds} --nimg=${winds} -k ${aforcecons};

    # Copy the newconstr.dat back to Xsede
    scp $pathn/optstring_calcs$step/newconstr.dat $xsedepath/s$step/
    scp $pathn/optstring_calcs$step/newconstr.RST $xsedepath/s$step/
    cd $pathn
fi


if [ "$2" == "copy" ]
then
    if [ -d $pathn/s$step.xsede ]
    then
        rm -r $pathn/s$step.xsede
    fi
    scp -r $xsedepath/tmp/ $pathn/s$step.xsede

    if [ -d $pathn/optstring_calcs$step ]
    then
        rm -r $pathn/optstring_calcs$step
    fi
    mkdir -p $pathn/optstring_calcs$step
fi

if [ "$2" == "caldis" ]
then
    #Atom numbers of each reaction coordiante are stored in the $pathn/dis_atnums.txt file
    for ((a=1;a<=${winds};a++));do
    echo "****Windows${a}****";
    cal_dis.py -i $pathn/string_atompairs.txt -p $pathn/s$step.xsede/${prmtop}.prmtop -c $pathn/s$step.xsede/f$a/${step}_production.netcdf -o $pathn/optstring_calcs$step/$a.dat -t $pathn/optstring_calcs$step/s${step}_i${a}_collect.dat;
    done
fi

if [ "$2" == "globin" ]
then
    #Global binning
    for ((a=1;a<=${dim};a++))
    do
        glo_bin.py -c $1 -w ${winds} -b 20 -d ${a}
    done
fi

if [ "$2" == "restr" ]
then
    # Regenerate the string
    cd $pathn/optstring_calcs$step;
    rege_string_test.py -i $pathn/string_atompairs.txt -d ${dim} --img=${winds} --nimg=${winds} -k ${aforcecons};
fi
   

if [ "$2" == "scopy" ]
then
    # Copy the newconstr.dat back to Xsede
    scp $pathn/optstring_calcs$step/newconstr.dat $xsedepath/s$step/
    scp $pathn/optstring_calcs$step/newconstr.RST $xsedepath/s$step/
    cd $pathn
fi


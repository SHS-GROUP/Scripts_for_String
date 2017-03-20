# script to prepare and submit jobs for QM/MM free energy simulations
############################################
# set STEP number of simulations
STEP=5;
PSTEP=$(($STEP-1));
#############################################
# set section
section="$1";
##################################################
# set pathway
IPATH="/oasis/projects/nsf/uic317/pfliuiuc/Projs/SLO/QMMM/test/string_running";
# In the ADIR (Charmm directory) we have sub directories for Charmm input files, slm files, and stream files
ADIR="/oasis/projects/nsf/uic317/pfliuiuc/Projs/SLO/QMMM/test/string_running/amber";
# set image number
num_img=20;
max_img=$(($num_img+1));
dim=3;
top="pre_react";
####################################################
# section to run the mm_equ step
if [ "$section" == "mmeq" ]
    then
    for(( i=1;i<max_img;i++)); do
        mkdir -p $IPATH/s$PSTEP/f$i/mm_equ;
        cd $IPATH/s$PSTEP/f$i/mm_equ/;
        cp $IPATH/interpol/f$i/${top}.prmtop .;
        cp $IPATH/interpol/f$i/interpol_${i}.inpcrd .;
	cp $ADIR/input_files/mm_equ.in .;
        cp $ADIR/slm_files/amber_mm_equ_xsede.slm mm_equ_xsede.slm;
        sed -i 's/IIIII/'$i'/g' mm_equ_xsede.slm;
	sed -i 's/TOPTOP/'$top'/g' mm_equ_xsede.slm;
        sbatch mm_equ_xsede.slm;
	sleep 5;
    done
fi
####################################################
# section to check whether the mm_equ finished
if [ "$section" == "mmeqcheck" ]
    then
    echo "This is a normal exit check" > $IPATH/s$PSTEP/mm_equ_finish.log;
    for(( i=1;i<max_img;i++)); do
	echo "***********IMAGE==$i*************" >> $IPATH/s$PSTEP/mm_equ_finish.log;
	tail -n 20 $IPATH/s$PSTEP/f$i/mm_equ/mm_equ.out >> $IPATH/s$PSTEP/mm_equ_finish.log;
    done
    echo "Run Done :"
    grep -i 'Run   done at' $IPATH/s$PSTEP/mm_equ_finish.log | wc -l
    echo "Wallclock Done :"
    grep -i 'wallclock' $IPATH/s$PSTEP/mm_equ_finish.log | wc -l
fi
####################################################
# section to run the qmmm_equ
if [ "$section" == "qmmmeq" ]
    then
    for(( i=1;i<max_img;i++)); do
        rm -r $IPATH/s$PSTEP/f$i/qmmm_equ;
        mkdir -p $IPATH/s$PSTEP/f$i/qmmm_equ;
        cd $IPATH/s$PSTEP/f$i/qmmm_equ/;
        cp $IPATH/interpol/f$i/${top}.prmtop .;
        cp $IPATH/s$PSTEP/f$i/mm_equ/mm_equ_f${i}.rst .;
        cp $ADIR/input_files/qmmm_equ.in .;
        cp $ADIR/input_files/qc_job.tpl .;
	sed -i 's/IIIII/'$i'/g' qmmm_equ.in;
        cp $IPATH/interpol/iniconstr.RST .;
        j=`echo "${i}*${dim}" | bc`
        head -$j iniconstr.RST | tail -${dim} > ${i}.RST;
        cp $ADIR/slm_files/amber_qmmm_equ_xsede.slm qmmm_equ_xsede.slm;
        sed -i 's/IIIII/'$i'/g' qmmm_equ_xsede.slm;
	sed -i 's/TOPTOP/'$top'/g' qmmm_equ_xsede.slm;
        sbatch qmmm_equ_xsede.slm;
	sleep 5;
    done
fi
############################################################
# section to check the completion of the qmmm_equ
# prepare the files for string update for the first iteration
# After it needs to go back to shscluster to run string_shs.sh
if [ "$section" == "qmmmeqcheck" ]
    then
    # Prepare the files required to update the string
    rm -r $IPATH/tmp;
    mkdir -p $IPATH/tmp;
    cp $IPATH/s0/f1/mm_equ/${top}.prmtop $IPATH/s$PSTEP/;
    cp $IPATH/s0/${top}.prmtop $IPATH/tmp/;
    echo "This is a normal exit check" > $IPATH/s$PSTEP/qmmm_equ_finish.log;
    for(( i=1;i<max_img;i++));do
        # Prepare files for string opitmization    
        mkdir -p $IPATH/tmp/f$i;
        cp $IPATH/s$PSTEP/f$i/qmmm_equ/qmmm_equ_f${i}.netcdf $IPATH/tmp/f$i/0_production.netcdf;
        # Generate a normal exit log file
        echo "***********IMAGE==$i*************" >> $IPATH/s$PSTEP/qmmm_equ_finish.log;
        tail -n 20 $IPATH/s$PSTEP/f$i/qmmm_equ/qmmm_equ.out >> $IPATH/s$PSTEP/qmmm_equ_finish.log;
    done
    echo "Run Done :"
    grep -i 'Run   done at' $IPATH/s$PSTEP/qmmm_equ_finish.log | wc -l
    echo "Wallclock Done :"
    grep -i 'wallclock' $IPATH/s$PSTEP/qmmm_equ_finish.log | wc -l
fi
#######################################################
# section to setup files for the first iteration
if [ "$section" == "submit1" ]
    then
    rm -rf $IPATH/s$STEP;
    for(( i=1;i<max_img;i++));do
        mkdir -p $IPATH/s$STEP/f$i;
        cd $IPATH/s$STEP/f$i;
        cp $IPATH/s0/f$i/qmmm_equ/${top}.prmtop .;
        cp $IPATH/s0/f$i/qmmm_equ/qmmm_equ_f${i}.rst ${STEP}_initial_f${i}.rst;
        cp $IPATH/s0/newconstr.RST newconstr.${STEP}.RST;
        echo $STEP > step;
	j=`echo "${i}*${dim}" | bc`
	head -$j newconstr.${STEP}.RST | tail -${dim} > ${i}.RST;
        cp $ADIR/input_files/qmmm_string.in .;
	sed -i 's/IIIII/'$i'/g' qmmm_string.in;
        cp $ADIR/input_files/qc_job.tpl .;
        cp $ADIR/slm_files/amber_qmmm_string_xsede.slm qmmm_string_xsede.slm;
	sed -i 's/SSSSS/'$STEP'/g' qmmm_string_xsede.slm;
	sed -i 's/IIIII/'$i'/g' qmmm_string_xsede.slm;
	sed -i 's/TOPTOP/'$top'/g' qmmm_string_xsede.slm;
        sbatch qmmm_string_xsede.slm
        sleep 5;
    done
    cd $IPATH/s$STEP;
    echo  " Step $STEP " >> ../runlog
    cat ./f1/newconstr.$STEP.RST >> ../runlog
fi
#######################################################
# After the 
# section to check whether the jobs are well finished
# and prepare the tmp directory sent to shsc cluster
if [ "$section" == "check" ]
    then
    rm -r $IPATH/tmp;
    mkdir -p $IPATH/tmp;
    cp $IPATH/s0/f1/mm_equ/${top}.prmtop $IPATH/s$PSTEP/;
    cp $IPATH/s0/${top}.prmtop $IPATH/tmp/;
    echo "This is a normal exit check" > $IPATH/s$PSTEP/qmmm_sampling_finish.log;
    for(( i=1;i<max_img;i++));do
        mkdir -p $IPATH/tmp/f$i;
        cp $IPATH/s$PSTEP/f$i/qmmm_string_${PSTEP}_f${i}.netcdf $IPATH/tmp/f$i/${PSTEP}_production.netcdf;
        echo "***********IMAGE==$i*************" >> $IPATH/s$PSTEP/qmmm_sampling_finish.log;
        tail -n 20 $IPATH/s$PSTEP/f$i/qmmm_string.out >> $IPATH/s$PSTEP/qmmm_sampling_finish.log;
    done
    echo "Run Done :"
    grep -i 'Run   done at' $IPATH/s$PSTEP/qmmm_sampling_finish.log | wc -l
    echo "Wallclock Done :"
    grep -i 'wallclock' $IPATH/s$PSTEP/qmmm_sampling_finish.log | wc -l
fi
###################################################################
# section to submit jobs
if [ "$section" == "submit" ]
    then
    for(( i=1;i<max_img;i++));do
        mkdir -p $IPATH/s$STEP/f$i;
        cd $IPATH/s$STEP/f$i;
        cp $IPATH/s$PSTEP/f$i/${top}.prmtop .;
        cp $IPATH/s$PSTEP/f$i/${PSTEP}_last_f${i}.rst ${STEP}_initial_f${i}.rst;
        cp $IPATH/s$PSTEP/newconstr.RST newconstr.${STEP}.RST;
        echo $STEP > step;
	j=`echo "${i}*${dim}" | bc`
	head -$j newconstr.${STEP}.RST | tail -${dim} > ${i}.RST;
        cp $ADIR/input_files/qmmm_string.in .;
        sed -i 's/IIIII/'$i'/g' qmmm_string.in;
        cp $ADIR/input_files/qc_job.tpl .;
	cp $ADIR/slm_files/amber_qmmm_string_xsede.slm qmmm_string_xsede.slm;
        sed -i 's/SSSSS/'$STEP'/g' qmmm_string_xsede.slm;
        sed -i 's/IIIII/'$i'/g' qmmm_string_xsede.slm;
	sed -i 's/TOPTOP/'$top'/g' qmmm_string_xsede.slm;
        sbatch qmmm_string_xsede.slm
        sleep 5;
    done
    cd $IPATH/s$STEP;
    echo  " Step $STEP " >> ../runlog
    cat ./f1/newconstr.$STEP.RST >> ../runlog
fi

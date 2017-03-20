# script to prepare and submit jobs for QM/MM free energy simulations
############################################
# set STEP number of simulations
STEP=3;
PSTEP=$(($STEP-1));
#############################################
# set section
section="$1";
##################################################
# set pathway
IPATH="/oasis/projects/nsf/uic317/pfliuiuc/Projs/String_test/4-mastro/transfer_to_maestro/running_test";
# In the CDIR (Charmm directory) we have sub directories for Charmm input files, slm files, and stream files
CDIR="/oasis/projects/nsf/uic317/pfliuiuc/Projs/String_test/4-mastro/transfer_to_maestro/running_test/charmm";
# set image number
num_img=18;
max_img=$(($num_img+1));
####################################################
# section to run the mm_equ step
if [ "$section" == "mmeq" ]
    then
    for(( i=1;i<max_img;i++)); do
        mkdir -p $IPATH/s$PSTEP/f$i/mm_equ;
        cd $IPATH/s$PSTEP/f$i/mm_equ/;
        cp -v $IPATH/interpol/f$i/react_link.psf .;
        cp -v $IPATH/interpol/f$i/interpol_${i}.crd .;
        cp -v $CDIR/input_files/qchem.in .;
	cp -v $CDIR/input_files/mm_equil.in .;
	sed -i 's/IIIII/'$i'/g' mm_equil.in;
        cp -v $CDIR/slm_files/charmm_mm_xsede.slm mm_equ_xsede.slm;
        sed -i 's/inputfile/mm_equil/g' mm_equ_xsede.slm;
        sed -i 's/IIIII/'$i'/g' mm_equ_xsede.slm;
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
	tail -n 20 $IPATH/s$PSTEP/f$i/mm_equ/mm_equil_f${i}.out >> $IPATH/s$PSTEP/mm_equ_finish.log;
	awk '{print $0}' $IPATH/s$PSTEP/f$i/mm_equ/mm_equil_f${i}.out | grep "DYNA>    10000" >> $IPATH/s$PSTEP/mm_equ_finish.log;
	awk '{print $0}' $IPATH/s$PSTEP/f$i/mm_equ/mm_equil_f${i}.out | grep "AVER>    10000" >> $IPATH/s$PSTEP/mm_equ_finish.log;
    done
fi
####################################################
# section to run the qmmm_equ
if [ "$section" == "qmmmeq" ]
    then
    for(( i=1;i<max_img;i++)); do
        mkdir -p $IPATH/s$PSTEP/f$i/qmmm_equ;
        cd $IPATH/s$PSTEP/f$i/qmmm_equ/;
        cp -v $IPATH/interpol/f$i/react_link.psf .;
        cp -v $IPATH/s$PSTEP/f$i/mm_equ/mm_equil.crd .;
        cp -v $CDIR/input_files/qchem.in .;
        cp -v $CDIR/input_files/qmmm_equil.in .;
	cp -v $IPATH/interpol/iniconstr.dat .;
	cp -v $IPATH/interpol/force.dat .;
        head -$i iniconstr.dat | tail -1 | awk '{ for ( i=1;i<=NF;i++ ) printf "SET V%s %6.4f \n",i,$i }' > val.stream;
        head -$i force.dat | tail -1 | awk '{ for ( i=1;i<=NF;i++ ) printf "SET F%s %6.4f \n",i,$i }' >> val.stream;
	cp -v $CDIR/stream_files/constr.stream .;
	cp -v $CDIR/stream_files/getcoord.stream .;
        cp -v $CDIR/slm_files/charmm_qmmm_xsede.slm qmmm_equ_xsede.slm;
        sed -i 's/inputfile/qmmm_equil/g' qmmm_equ_xsede.slm;
        sed -i 's/IIIII/'$i'/g' qmmm_equ_xsede.slm;
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
    cp -v $IPATH/s0/f1/mm_equ/react_link.psf $IPATH/s$PSTEP/;
    cp -v $IPATH/s0/react_link.psf $IPATH/tmp/;
    echo "This is a normal exit check" > $IPATH/s$PSTEP/qmmm_equ_finish.log;
    for(( i=1;i<max_img;i++));do
        # Prepare files for string opitmization    
        mkdir -p $IPATH/tmp/f$i;
        cp $IPATH/s$PSTEP/f$i/qmmm_equ/qmmm_equil.dcd $IPATH/tmp/f$i/0_production.dcd;
        # Generate a normal exit log file
        echo "***********IMAGE==$i*************" >> $IPATH/s$PSTEP/qmmm_equ_finish.log;
        tail -n 20 $IPATH/s$PSTEP/f$i/qmmm_equ/qmmm_equil_f${i}.out >> $IPATH/s$PSTEP/qmmm_equ_finish.log;
        awk '{print $0}' $IPATH/s$PSTEP/f$i/qmmm_equ/qmmm_equil_f${i}.out | grep "DYNA>      100" >> $IPATH/s$PSTEP/qmmm_equ_finish.log;
        awk '{print $0}' $IPATH/s$PSTEP/f$i/qmmm_equ/qmmm_equil_f${i}.out | grep "AVER>      100" >> $IPATH/s$PSTEP/qmmm_equ_finish.log;
    done
fi
#########################################################
# section to get last coordinates of the qmmm_equ
if [ "$section" == "getfinalcrd0" ]
    then
    for(( i=1;i<max_img;i++));do
        cd $IPATH/s$PSTEP/f$i/qmmm_equ;
        cp -v $CDIR/input_files/get_final_coord.in .;
	sed -i 's/initialcoord/qmmm_equil/g' get_final_coord.in;
        sed -i 's/trajectory.dcd/qmmm_equil.dcd/g' get_final_coord.in;
	sed -i 's/SSSSS/'${PSTEP}'/g' get_final_coord.in;
        cp -v $CDIR/slm_files/charmm_mm_xsede.slm get_final_coord_xsede.slm;
        sed -i 's/inputfile/get_final_coord/g' get_final_coord_xsede.slm;
        sed -i 's/IIIII/'$i'/g' get_final_coord_xsede.slm;
        sbatch get_final_coord_xsede.slm;
        sleep 5;
    done
fi
#######################################################
# section to setup files for the first iteration
if [ "$section" == "submit1" ]
    then
    for(( i=1;i<max_img;i++));do
        mkdir -p $IPATH/s$STEP/f$i;
        cd $IPATH/s$STEP/f$i;
        cp -v $IPATH/s0/f$i/qmmm_equ/react_link.psf .;
        cp -v $IPATH/s0/f$i/qmmm_equ/${PSTEP}_last.crd ${STEP}_initial.crd;
        cp -v $IPATH/s0/newconstr.dat newconstr.$STEP.dat;
        cp -v $IPATH/interpol/force.dat force.$STEP.dat;
        echo $STEP > step;
        head -$i newconstr.$STEP.dat | tail -1 | awk '{ for ( i=1;i<=NF;i++ ) printf "SET V%s %6.4f \n",i,$i }' > val.stream;
        head -$i force.$STEP.dat | tail -1 | awk '{ for ( i=1;i<=NF;i++ ) printf "SET F%s %6.4f \n",i,$i }' >> val.stream;
        cp -v $CDIR/input_files/qchem.in .;
        cp -v $CDIR/stream_files/*.stream .;
        cp -v $CDIR/input_files/qmmm_sampling.in .;
        sed -i 's/SSSSS/'${STEP}'/g' qmmm_sampling.in;
        cp -v $CDIR/slm_files/charmm_qmmm_xsede.slm qmmm_sampling_xsede.slm;
        sed -i 's/inputfile/qmmm_sampling/g' qmmm_sampling_xsede.slm;
	sed -i 's/IIIII/'$i'/g' qmmm_sampling_xsede.slm;
        sbatch qmmm_sampling_xsede.slm;
        sleep 5;
    done
    cd $IPATH/s$STEP;
    echo  " Step $STEP " >> ../runlog
    cp -v $IPATH/s0/newconstr.dat newconstr.$STEP.dat;
    cp -v $IPATH/interpol/force.dat force.$STEP.dat;
    cat newconstr.$STEP.dat >> ../runlog
    cat force.$STEP.dat >> ../runlog
fi
#######################################################
# After the 
# section to check whether the jobs are well finished
# and prepare the tmp directory sent to shsc cluster
if [ "$section" == "check" ]
    then
    rm -r $IPATH/tmp;
    mkdir -p $IPATH/tmp;
    cp -v $IPATH/s0/f1/mm_equ/react_link.psf $IPATH/s$PSTEP/;
    cp -v $IPATH/s0/react_link.psf $IPATH/tmp/;
    echo "This is a normal exit check" > $IPATH/s$PSTEP/qmmm_sampling_finish.log;
    for(( i=1;i<max_img;i++));do
        mkdir -p $IPATH/tmp/f$i;
        cp $IPATH/s$PSTEP/f$i/${PSTEP}_production.dcd $IPATH/tmp/f$i/;
        echo "***********IMAGE==$i*************" >> $IPATH/s$PSTEP/qmmm_sampling_finish.log;
        tail -n 5 $IPATH/s$PSTEP/f$i/qmmm_sampling_f${i}.out >> $IPATH/s$PSTEP/qmmm_sampling_finish.log;
        awk '{print $0}' $IPATH/s$PSTEP/f$i/qmmm_sampling_f${i}.out | grep "DYNA>      100" >> $IPATH/s$PSTEP/qmmm_sampling_finish.log;
        awk '{print $0}' $IPATH/s$PSTEP/f$i/qmmm_sampling_f${i}.out | grep "AVER>      100" >> $IPATH/s$PSTEP/qmmm_sampling_finish.log;
    done
fi
###########################################################
# section to get last coordinates of the previous batch of jobs
if [ "$section" == "getfinalcrd" ]
    then
    for(( i=1;i<max_img;i++));do
        cd $IPATH/s$PSTEP/f$i/;
	cp -v $CDIR/input_files/get_final_coord.in .;
	sed -i 's/initialcoord.crd/step'${PSTEP}'.cor/g' get_final_coord.in;
        sed -i 's/trajectory.dcd/'${PSTEP}'_production.dcd/g' get_final_coord.in;
	sed -i 's/SSSSS/'${PSTEP}'/g' get_final_coord.in;
        cp -v $CDIR/slm_files/charmm_mm_xsede.slm get_final_coord_xsede.slm;
	sed -i 's/inputfile/get_final_coord/g' get_final_coord_xsede.slm;
	sed -i 's/IIIII/'$i'/g' get_final_coord_xsede.slm;
	sbatch get_final_coord_xsede.slm;
        sleep 5;
    done
fi
###################################################################
# section to submit jobs
if [ "$section" == "submit" ]
    then
    for(( i=1;i<max_img;i++));do
        mkdir -p $IPATH/s$STEP/f$i;
        cd $IPATH/s$STEP/f$i;
        cp -v $IPATH/s$PSTEP/f$i/react_link.psf .;
        cp -v $IPATH/s$PSTEP/f$i/${PSTEP}_last.crd ${STEP}_initial.crd;
        cp -v $IPATH/s$PSTEP/newconstr.dat newconstr.$STEP.dat;
        cp -v $IPATH/interpol/force.dat force.$STEP.dat;
        echo $STEP > step;
        head -$i newconstr.$STEP.dat | tail -1 | awk '{ for ( i=1;i<=NF;i++ ) printf "SET V%s %6.4f \n",i,$i }' > val.stream;
        head -$i force.$STEP.dat | tail -1 | awk '{ for ( i=1;i<=NF;i++ ) printf "SET F%s %6.4f \n",i,$i }' >> val.stream;
        cp -v $CDIR/input_files/qchem.in .;
        cp -v $CDIR/stream_files/*.stream .;
        cp -v $CDIR/input_files/qmmm_sampling.in .;
        sed -i 's/SSSSS/'${STEP}'/g' qmmm_sampling.in;
        cp -v $CDIR/slm_files/charmm_qmmm_xsede.slm qmmm_sampling_xsede.slm;
        sed -i 's/inputfile/qmmm_sampling/g' qmmm_sampling_xsede.slm;
        sed -i 's/IIIII/'$i'/g' qmmm_sampling_xsede.slm;
        sbatch qmmm_sampling_xsede.slm;
        sleep 5;
        done
    cd $IPATH/s$STEP;
    echo  " Step $STEP " >> ../runlog
    cp -v $IPATH/s$PSTEP/newconstr.dat newconstr.$STEP.dat;
    cp -v $IPATH/interpol/force.dat force.$STEP.dat;
    cat newconstr.$STEP.dat >> ../runlog
    cat force.$STEP.dat >> ../runlog
fi

#! /bin/bash

NP=3

make rectangular_plate_fvk

# Stuff to move
important_files="rectangular_plate_foeppl_von_karman.cc
                 rectangular_plate_fvk
                 rectangular_run.bash   "

# Setup directories YOU MUST PICK A NAME FOR YOUR OURPUT DIRECTORY.
main_dir=NEW_TEST_LU_rectangular_np_${NP}
if [ -e $main_dir ]; then
    echo " "
    echo "WARNING: Directory " $main_dir "already exists!"
    read -p "         remove it and continue? [Y/n] " yn
    case $yn in
        ''|[Yy]* ) rm -rf $main_dir; break;;
        [Nn]* ) echo "Can't continue until you move $main_dir"; exit;;
    esac
    echo " "
fi
mkdir $main_dir

cp $important_files $main_dir
cd $main_dir



for i in {1..1}
do
	#ETA=1.0e6
	NU=0.495
	MU=1.0e6

	reslt_dir=RESLT
	mkdir $reslt_dir

	#-------------------------------------------------------------
	# Oomph-lib time!
	#-------------------------------------------------------------
	echo "Doing "$reslt_dir
	mpirun -np $NP ./rectangular_plate_fvk \
	  --dir $reslt_dir  \
	  --mu  $MU         \
          2>&1 | tee $reslt_dir/OUTPUT
done

echo " "
echo " "
echo "Done!"
echo " "
echo " "


exit

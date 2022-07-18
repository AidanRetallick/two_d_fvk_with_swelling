#! /bin/bash

make axisym_displacement_control_fvk_microphone_with_contact

# Stuff to move
important_files="axisym_displacement_control_fvk_microphone_with_contact.cc
                 axisym_displacement_control_fvk_microphone_with_contact "

# Setup directories YOU MUST PICK A NAME FOR YOUR OURPUT DIRECTORY.
main_dir=0   # Change 0 to be the name of your output directory. e.g. NEW_xyz
if [ main_dir == 0 ]
   echo " "
   echo "You need to give your new output a name"
fi    
if [ -e $main_dir ]; then
    echo " "
    echo "WARNING: Directory " $main_dir "already exists -- removing it."
    echo " "
    rm -rf $main_dir
fi
mkdir $main_dir

cp $important_files $main_dir
cd $main_dir

reslt_dir=RESLT
mkdir $reslt_dir

#-------------------------------------------------------------
# Oomph-lib time!
#-------------------------------------------------------------
echo "Doing "$reslt_dir
./axisym_displacement_control_fvk_microphone_with_contact \
                          --dir $reslt_dir \
                              > $reslt_dir/OUTPUT

echo " "
echo " "
echo "Done!"
echo " "
echo " "


exit

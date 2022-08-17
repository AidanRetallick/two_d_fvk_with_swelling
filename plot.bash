#!/bin/bash

echo " "

SAVE=false
TRUESCALE=false
TRUESCALE_STRING=""
X=1
Y=2
XLABEL="Time"
YLABEL="Centre deflection, w_0"

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	-d|--directory)
	    echo "  $key was passed, expecting NEW_ directory" 
	    shift
            echo "  plotting from the RESLT_ directories in ${1}"
            DIRECTORY=$1
            shift
	    ;;
	-s|--save)
	    echo "  $key was passed, saving file."
	    SAVE=true
            shift
            SAVENAME=$1
	    shift
	    ;;
	*)
	    echo "ERROR!  $key is an invalid argument!"
	    echo " "
	    exit
    esac
done

FILENAMEARRAY=$(echo ${DIRECTORY}/RESLT*/trace_nondim.dat)
echo $FILENAMEARRAY
gnuplot -p << EOF
FILES="${FILENAMEARRAY[*]}"

if ("$SAVE" eq "true") {
 print " "
 print "Saving Plot as $SAVENAME"
 set terminal postscript eps size 5,2.6 enhanced color font 'Helvetica,16'
 set output "Plots/${SAVENAME}"
}
 set xlabel "$XLABEL"
 set ylabel "$YLABEL"
 set key on right bottom
 set xrange [0:100]
 set yrange [0:0.16]

  plot for [i=1:words(FILES)] word(FILES,i) using 4:3 with linespoints lt i title "mu = 1.0e4, eta = 1.0e".i, \
       for [i=1:words(FILES)] word(FILES,i) using 4:3:5 every 2 with labels tc lt i offset 0.5,-0.7 nopoint notitle

pause -1
EOF

echo " "

exit

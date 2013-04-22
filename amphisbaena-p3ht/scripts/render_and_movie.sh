for i in *.pov
do
 png=` echo $i | sed s/pov/png/ `

 if [ ! -f "${png}" ]
 then
       	povray -D0 -W800 -H800 "${i}"
 else
	 echo "${i} : ${png} Appears to be rendered... skipping"
 fi

done

mencoder "mf://*.png" -ovc lavc -lavcopts vcodec=mjpeg -fps 12 -o pretty.avi


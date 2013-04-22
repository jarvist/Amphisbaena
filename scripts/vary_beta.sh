#for  in 0.00 0.001 0.01 0.05 0.1
#do
for dense in 0.175 0.20 0.225 0.25 0.275 0.30  0.35 0.40 0.45 0.50
do
 echo $bias $dense
 
  cat amphisbaena.h.bak |  sed s/dense/${dense}/ | sed s/mountain/50/ > amphisbaena.h 
#  cat gorgophone.c.bak | sed s/bias/${bias}/ > gorgophone.c
  gcc -o a_${dense}.bin amphisbaena.c -lm -O4
  
  cat amphisbaena.h.bak | sed s/dense/${dense}/ | sed s/mountain/-50/ > amphisbaena.h
 # cat gorgophone.c.bak | sed s/bias/${bias}/ > gorgophone.c  
  gcc -o r_${dense}.bin amphisbaena.c -lm -O4
      
done
#done


for dense in 0.60 0.70 0.80
#0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50
#0.26 0.27 0.28 0.29 0.31 0.32 0.33 0.34
#0.16 0.17 0.18 0.19 0.21 0.22 0.23 0.24 0.25
#0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50
do
 echo $dense
 
  cat amphisbaena.h.bak |  sed s/dense/${dense}/ | sed s/mountain/50/ > amphisbaena.h 
#  cat gorgophone.c.bak | sed s/bias/${bias}/ > gorgophone.c
  gcc -o a_${dense}.bin amphisbaena.c -lm -O4
  
  cat amphisbaena.h.bak | sed s/dense/${dense}/ | sed s/mountain/-50/ > amphisbaena.h
 # cat gorgophone.c.bak | sed s/bias/${bias}/ > gorgophone.c  
  gcc -o r_${dense}.bin amphisbaena.c -lm -O4
      
done



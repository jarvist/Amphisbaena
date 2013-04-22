read bias

for bigrock in  0.6 0.45 0.40 0.35 0.30 0.25 0.20 0.575 0.525 0.475 0.425 0.375 0.325 0.275 0.225
do
 cat amphisbaena.h.bak | sed s/bigrock/${bigrock}/ | sed s/candy/50/ > amphisbaena.h
 cat gorgophone.c.bak | sed s/bias/${bias}/ > gorgophone.c
 gcc -o attractive.bin amphisbaena.c -lm -O4
 sh ./attractive.bin > a_${bias}_${bigrock}.dat
 
 cat amphisbaena.h.bak | sed s/bigrock/${bigrock}/ | sed s/candy/-50/ > amphisbaena.h
 cat gorgophone.c.bak | sed s/bias/${bias}/ > gorgophone.c
 gcc -o repulsive.bin amphisbaena.c -lm -O4
 ./repulsive.bin > r_${bias}_${bigrock}.dat

done

for i in snakes_*.df3
do
 cp "${i}" snakes.df3
 povray -W600 -H600 "-O${i}.png" ToF_density.pov
 mv "${i}" old/
 read p
done

max=171

for i in *.df3
do
 num=` echo $i | sed -e s/snakes_// -e s/\.df3// `
 cat ToF_density.pov | sed -e "s/clock/(${num}\/${max})/g" -e "s/snakes.df3/${i}/" > ${i}.pov 
done

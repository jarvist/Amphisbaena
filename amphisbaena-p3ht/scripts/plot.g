set logscale
set output "pcbm_vs_p3ht.eps"
set terminal postscript enhanced color
p "repulse/end/pcbm/x.tof" u 2:3 w lp, "repulse/end/p3ht/x.tof" u 2:3 w lp, "aggregate/1000MS/pcbm/x.tof" u 2:3 w lp, "aggregate/1000MS/p3ht/x.tof" u 2:3 w lp, "Ea=Eb=Ec/S00000032/pcbm/x.tof" u 2:3 w lp, "Ea=Eb=Ec/S00000032/p3ht/x.tof" u 2:3 w lp

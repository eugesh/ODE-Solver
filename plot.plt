# splot "rk.dat" using 2:3:4 w l, "ae.dat" using 2:3:4 w l, "ai.dat" using 2:3:4 w l
plot "rosen.dat" using 2:3 w l, "rk.dat" using 2:3 w l, "ae.dat" using 2:3 w l, "ai.dat" using 2:3 w l, "precor.dat" using 2:3 w l
pause -1

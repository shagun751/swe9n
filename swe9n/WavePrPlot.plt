set terminal pngcairo size 1280, 640 linewidth 3
#set yrange [0.85:1.15]
#set xrange [0:45]
#set ytics 0.96,0.005,1.04
#set ytics 0.02
set key left top
set xtics 2
set grid ytics lt 0 lw 1
set grid xtics lt 0 lw 1
show grid
set palette model HSV
do for [ii=1:2] {
	fileN3=sprintf("WaveProbe_v5.34.tnF.tB.w2.dat")
	fileN2=sprintf("WaveProbe_v5.34e.tnF.tB.dat")
	fileN1=sprintf("WaveProbe_tnF.dat")
	fileN=sprintf("WEta%d.png",ii)	
	
	set output fileN
	iib=1+(ii-1)*4+2
	plot fileN3 using ($1/3600):iib with lines, fileN2 using ($1/3600):iib with lines, fileN1 using ($1/3600):iib with lines	

	fileN=sprintf("WP%d.png",ii)	
	set output fileN
	iib=1+(ii-1)*4+3
	plot fileN3 using ($1/3600):iib with lines, fileN2 using ($1/3600):iib with lines, fileN1 using ($1/3600):iib with lines	

	fileN=sprintf("WQ%d.png",ii)	
	set output fileN
	iib=1+(ii-1)*4+3
	plot fileN3 using ($1/3600):iib with lines, fileN2 using ($1/3600):iib with lines, fileN1 using ($1/3600):iib with lines	

}

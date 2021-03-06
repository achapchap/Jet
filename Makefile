#-----------------------------------     
#    Makefile per SPLASH 
#-----------------------------------     
F77=gfortran 
############## Linux
SWI= -c -O3 -Wall -fcheck=bounds -llapack -lblas
COM= -O3 -Wall -fcheck=bounds -llapack -lblas

moduli=	main.o\
	bcpre.o\
	calsol.o\
	calsolp.o\
	calvel.o\
	check.o\
	delta.o\
	distri.o\
	disun2.o\
	doldfil1.o\
	filiniz.o\
	finte.o\
	gauss.o\
	get.o\
	initial.o\
	input.o\
	intret.o\
	nortan.o\
	poli3.o\
	prefin.o\
	ridis6.o\
	ropn.o\
	rot.o\
	shallo.o\
	solv1.o\
	solv2.o\
	solv22.o\
	solv22p.o\
	solv2p.o\
	spline.o\
	splint.o\
	spline1.o\
	splint1.o\
	spline3.o\
	splint3.o\
	splver2.o\
	stampa.o

spl: $(moduli)
	$(F77) $(moduli) $(COM)   -o 2dt

main.o     : main.f90
	$(F77) $(SWI) main.f90
bcpre.o     : bcpre.f
	$(F77) $(SWI) bcpre.f
calsol.o     : calsol.f
	$(F77) $(SWI) calsol.f
calsolp.o     : calsolp.f
	$(F77) $(SWI) calsolp.f
calvel.o     : calvel.f
	$(F77) $(SWI) calvel.f
check.o     : check.f
	$(F77) $(SWI) check.f
delta.o     : delta.f
	$(F77) $(SWI) delta.f
distri.o     : distri.f
	$(F77) $(SWI) distri.f
disun2.o     : disun2.f
	$(F77) $(SWI) disun2.f
doldfil1.o     : doldfil1.f
	$(F77) $(SWI) doldfil1.f
filiniz.o     : filiniz.f
	$(F77) $(SWI) filiniz.f
finte.o     : finte.f
	$(F77) $(SWI) finte.f
gauss.o     : gauss.f
	$(F77) $(SWI) gauss.f
get.o     : get.f
	$(F77) $(SWI) get.f
initial.o     : initial.f
	$(F77) $(SWI) initial.f
input.o     : input.f
	$(F77) $(SWI) input.f
intret.o     : intret.f
	$(F77) $(SWI) intret.f
nortan.o     : nortan.f
	$(F77) $(SWI) nortan.f
poli3.o     : poli3.f
	$(F77) $(SWI) poli3.f
prefin.o     : prefin.f
	$(F77) $(SWI) prefin.f
ridis6.o     : ridis6.f90
	$(F77) $(SWI) ridis6.f90
ropn.o     : ropn.f
	$(F77) $(SWI) ropn.f
rot.o     : rot.f
	$(F77) $(SWI) rot.f
shallo.o     : shallo.f90
	$(F77) $(SWI) shallo.f90
solv1.o     : solv1.f
	$(F77) $(SWI) solv1.f
solv2.o     : solv2.f
	$(F77) $(SWI) solv2.f
solv22.o     : solv22.f90
	$(F77) $(SWI) solv22.f90
solv22p.o     : solv22p.f
	$(F77) $(SWI) solv22p.f
solv2p.o     : solv2p.f
	$(F77) $(SWI) solv2p.f
spline.o     : spline.f
	$(F77) $(SWI) spline.f
splint.o     : splint.f
	$(F77) $(SWI) splint.f
spline1.o     : spline1.f
	$(F77) $(SWI) spline1.f
splint1.o     : splint1.f
	$(F77) $(SWI) splint1.f
spline3.o     : spline3.f
	$(F77) $(SWI) spline3.f
splint3.o     : splint3.f
	$(F77) $(SWI) splint3.f
splver2.o     : splver2.f
	$(F77) $(SWI) splver2.f
stampa.o     : stampa.f
	$(F77) $(SWI) stampa.f

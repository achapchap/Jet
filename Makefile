#-----------------------------------     
#    Makefile per SPLASH 
#-----------------------------------     
Cfor=gfortran 
############## Linux
# SWI= -c -O3 -Wall -fcheck=bounds -llapack -lblas
# COM= -O3 -Wall -fcheck=bounds -llapack -lblas

SWI= -c -O3 -fcheck=all -ffpe-trap=invalid,zero  -fbacktrace -C -g3 -fcheck=bounds -llapack -lblas 
COM= -O3 -fcheck=all -ffpe-trap=invalid,zero  -fbacktrace -C -g3 -fcheck=bounds -llapack -lblas

moduli=	parameters.o\
	main.o\
	input.o\
	initial.o\
	filiniz.o\
	spline.o\
	spline1.o\
        spline3.o\
        splint.o\
        splint3.o\
        delta.o\
	solver.o\
	finte.o\
	caveta.o\
	check.o\
        splver.o\
        intret.o\
        poli3.o\
	ridis6.o\
        nortan.o\
        disuni.o\
        shallo.o\
        doldfil1.o\
        output.o\
        dipole.o\
        get.o
         
      
spl: $(moduli)
	$(Cfor) $(COM) $(moduli) -o spl -llapack -lblas 

parameters.o     : parameters.f90
	$(Cfor) $(SWI) parameters.f90
main.o     : main.f90
	$(Cfor) $(SWI) main.f90
input.o     : input.f90
	$(Cfor) $(SWI) input.f90
initial.o     : initial.f90
	$(Cfor) $(SWI) initial.f90
filiniz.o     : filiniz.f90
	$(Cfor) $(SWI) filiniz.f90
spline.o      : spline.f90
	$(Cfor) $(SWI) spline.f90
spline1.o      : spline1.f90
	$(Cfor) $(SWI) spline1.f90
spline3.o      : spline3.f90
	$(Cfor) $(SWI) spline3.f90
splint.o      : splint.f90
	$(Cfor) $(SWI) splint.f90
splint3.o      : splint3.f90
	$(Cfor) $(SWI) splint3.f90
delta.o      : delta.f90
	$(Cfor) $(SWI) delta.f90
solver.o      : solver.f90
	$(Cfor) $(SWI) solver.f90
finte.o       : finte.f90
	$(Cfor) $(SWI) finte.f90
caveta.o      : caveta.f90
	$(Cfor) $(SWI) caveta.f90
check.o      : check.f90
	$(Cfor) $(SWI) check.f90
splver.o      : splver.f90
	$(Cfor) $(SWI) splver.f90
intret.o      : intret.f90
	$(Cfor) $(SWI) intret.f90
poli3.o       : poli3.f90
	$(Cfor) $(SWI) poli3.f90
ridis6.o       : ridis6.f90
	$(Cfor) $(SWI) ridis6.f90
nortan.o       : nortan.f90
	$(Cfor) $(SWI) nortan.f90
disuni.o       : disuni.f90
	$(Cfor) $(SWI) disuni.f90
shallo.o       : shallo.f90
	$(Cfor) $(SWI) shallo.f90
doldfil1.o       : doldfil1.f90
	$(Cfor) $(SWI) doldfil1.f90
output.o       : output.f90
	$(Cfor) $(SWI) output.f90
dipole.o       : dipole.f90
	$(Cfor) $(SWI) dipole.f90
get.o       : get.f90
	$(Cfor) $(SWI) get.f90

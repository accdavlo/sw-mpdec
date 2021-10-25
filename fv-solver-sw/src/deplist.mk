equation.o: equation.f90 \
	finitevolume_vars.o \
	reconstruction.o

finitevolume_vars.o: finitevolume_vars.f90

jacobi_iteration.o: jacobi_iteration.f90

finitevolume.o: finitevolume.f90 \
	equation.o \
	finitevolume_vars.o \
	reconstruction.o \
	shocksindicator.o \
	mesh.o

mesh.o: mesh.f90 \
	finitevolume_vars.o 

output.o: output.f90 \
	equation.o \
	finitevolume_vars.o 

parameters.o: parameters.f90 \
	finitevolume_vars.o 

reconstruction.o: reconstruction.f90 \
	finitevolume_vars.o 

shocksindicator.o: shocksindicator.f90 \
	finitevolume_vars.o 

timediscretization.o: timediscretization.f90 \
	equation.o \
	finitevolume.o \
	finitevolume_vars.o \
	output.o \
	jacobi_iteration.o

postprocessing.o: postprocessing.f90 \
	equation.o \
	finitevolume_vars.o 

main.o: main.f90 \
	finitevolume.o \
	mesh.o \
	parameters.o \
	postprocessing.o\
	timediscretization.o 

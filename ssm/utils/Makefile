program: nrtype.f90 nrutil.f90 random.f90
	 gfortran -c -fPIC likfort2count.f90 nrtype.f90 nrutil.f90 random.f90 
	 R CMD SHLIB likfort2count.o random.o nrutil.o nrtype.o -lblas -llapack


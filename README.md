# Multivariate models for count data

Basic implementation of state-space and log-linear models for count data (bivariate case). 

Cite: 
@article{shapovalova2021multivariate,
  title={Multivariate count data models for time series forecasting},
  author={Shapovalova, Yuliya and Ba{\c{s}}t{\"u}rk, Nalan and Eichler, Michael},
  journal={Entropy},
  volume={23},
  number={6},
  pages={718},
  year={2021},
  publisher={MDPI}
}

## Models

ssm -- directory with the state-space model example; 

log-linear -- directory with the log-linear model example;

Both models contain and example file called main.R in corresponding directories.

## Dependencies
SSM depends on Fortran routine for SMC; you likely will need to recompile likfort2count.f90 function. Go to ssm/utils/ and run make from Terminal/command line. Makefile may need to be adjusted depending on your system.

log-linear model depends on C function; you likely will need to recompile gpoisson.c. Go to log-linear/utils and run make from Terminal/command line. Makefile may need to be adjusted depending on your system.










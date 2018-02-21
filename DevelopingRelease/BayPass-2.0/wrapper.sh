#!/bin/bash

##tar -zxvf baypass_2.1.tar.gz

##make clean all FC=ifort

chmod +x sources/i_baypass
ARGS="${output}"

##General Options
if [ -z "${npop}" ]; then
	echo "number of populations not set"
else
	echo "number of populations set to ${npop}"
	ARGS="${ARGS} -npop ${npop}"
fi
if [ -z "${gfile}" ]; then
	echo "genotyping data file not set"
else
	echo "genotyping data file set to ${gfile}"
	ARGS="${ARGS} -gfile ${gfile}"
fi
if [ -z "${efile}" ]; then
	echo "covariate mode not activated"
else
	echo "covariate file set to ${efile}"
	ARGS="${ARGS} -efile ${efile}"
fi
if [ -z "${scalecov}" ]; then
	echo "scale covariates not set"
else 
	echo "scale covariates set to ${scalecov}"
	ARGS="${ARGS} -scalecov ${scalecov}"
fi
if [ -z ${poolsizefile} ]; then
	echo "poolseq mode not activated"
else
	echo "Pool Size file is set to ${poosizefile}"
	ARGS="${ARGS} -poolsizefile ${poolsizefile}"
fi
if [ -z ${outprefix} ]; then
	echo "prefix for output file not set"
else
	echo "prefix for output file is set to ${outprefix}"
	ARGS="${ARGS} -outprefix ${outprefix}"
fi
##Modeling Options
if [ -z ${omegafile} ]; then
	echo "omega matrix file file not set"
else
	echo "omega matrix file is set to ${omegafile}"
	ARGS="${ARGS} -omegafile ${omegafile}"
fi
if [ -z ${rho} ]; then
	echo "rho parameter not set"
else
	echo "Rho parameter is set to ${rho}"
	ARGS="${ARGS} -rho ${rho}"
fi
if [ -z ${setpibetapar} ]; then
	echo "Pi beta prior parameter estimation inactivated"
else
	echo "Estimation of Pi beta prior active"
	ARGS="${ARGS} -setpibetaprior ${setpibetaprior}"
fi
if [ -z "${betapipriormin}" ]; then
	echo "Pi beta prior minimum parameter not set"
else 
	echo "Pi beta prior minimum parameter set to ${betapipriormin}"
	ARGS="${ARGS} -betapipriormin ${betapipriormin}"
fi
if [ -z "${betapipriormax}" ]; then
	echo "Pi beta prior maximum parameter not set"
else 
	echo "Pi beta prior maximum parameter set to ${betapipriormax}"
	ARGS="${ARGS} -betapipriormax ${betapipriormax}"
fi
if [ -z "${minbeta}" ]; then
	echo "Lower beta coefficients not set"
else 
	echo "Lower beta coefficients set to ${minbeta}"
	ARGS="${ARGS} -minbeta ${minbeta}"
fi
if [ -z "${maxbeta}" ]; then
	echo "Upper beta coefficients not set"
else 
	echo "Upper beta coefficients set to ${maxbeta}"
	ARGS="${ARGS} -maxbeta ${maxbeta}"
fi
##IS Covariate mode
if [ -z "${nbetagrid}" ]; then
	echo "IS covariate mode grid points not set"
else 
	echo "IS covariate mode grid points set to ${nbetagrid}"
	ARGS="${ARGS} -nbetagrid ${nbetagrid}"
fi
##MCMC Covariate mode
if [ -z "${covmcmc}" ]; then
	echo "MCMC covariate mode not activated"
else 
	echo "MCMC covariate mode activated"
	ARGS="${ARGS} -covmcmc ${covmcmc}"
fi
if [ -z "${auxmodel}" ]; then
	echo "Auxiliary variable mode not activated"
else 
	echo "Auxiliary model mode activated"
	ARGS="${ARGS} -auxmodel ${auxmodel}"
fi
if [ -z "${isingbeta}" ]; then
	echo "Beta of the Ising model not set"
else 
	echo "Beta of the Ising model set to ${isingbeta}"
	ARGS="${ARGS} -isingbeta ${isingbeta}"
fi
if [ -z "${auxPbetapriormin}" ]; then
	echo "auxiliary P beta prior minimum parameter not set"
else 
	echo "auxiliary P beta prior minimum parameter set to ${auxPbetapriormin}"
	ARGS="${ARGS} -auxPbetapriormin ${auxPbetapriormin}"
fi
if [ -z "${auxPbetapriormax}" ]; then
	echo "auxiliary P beta prior maximum parameter not set"
else 
	echo "auxiliary P beta prior maximum parameter set to ${auxPbetapriormax}"
	ARGS="${ARGS} -auxPbetapriormax ${auxPbetapriormax}"
fi
##MCMC Options
if [ -z "${nthreads}" ]; then
	echo "Number of threads for MCMC not set"
else 
	echo "Number of threads for MCMC set to ${nthreads}"
	ARGS="${ARGS} -nthreads ${nthreads}"
fi
if [ -z "${nval}" ]; then
	echo "Number of post-burnin and thinned samples to generate for MCMC not set"
else 
	echo "Number of post-burnin and thinned samples to generate for MCMC set to ${nval}"
	ARGS="${ARGS} -nval ${nval}"
fi
if [ -z "${thin}" ]; then
	echo "Size of thining for MCMC not set"
else 
	echo "Size of thining for MCMC set to ${thin}"
	ARGS="${ARGS} -thin ${thin}"
fi
if [ -z "${burnin}" ]; then
	echo "Burn-in length for MCMC not set"
else 
	echo "Burn-in length for MCMC set to ${burnin}"
	ARGS="${ARGS} -burnin ${burnin}"
fi
if [ -z "${npilot}" ]; then
	echo "Number of pilot runs for MCMC not set"
else 
	echo "Number of pilot runs for MCMC set to ${npilot}"
	ARGS="${ARGS} -npilot ${npilot}"
fi
if [ -z "${pilotlength}" ]; then
	echo "Number of pilot runs for MCMC not set"
else 
	echo "Number of pilot runs for MCMC set to ${pilotlength}"
	ARGS="${ARGS} -pilotlength ${pilotlength}"
fi
if [ -z "${accinf}" ]; then
	echo "Lower target acceptance rate bound for MCMC not set"
else 
	echo "Lower target acceptance rate bound for MCMC set to ${accinf}"
	ARGS="${ARGS} -accinf ${accinf}"
fi
if [ -z "${accsup}" ]; then
	echo "Upper target acceptance rate bound for MCMC not set"
else 
	echo "Upper target acceptance rate bound for MCMC set to ${accsup}"
	ARGS="${ARGS} -accsup ${accsup}"
fi
if [ -z "${adjrate}" ]; then
	echo "Adjustment factor for MCMC not set"
else 
	echo "Adjustment factor for MCMC set to ${adjrate}"
	ARGS="${ARGS} -adjrate ${adjrate}"
fi
if [ -z "${adjrate}" ]; then
	echo "Adjustment factor for MCMC not set"
else 
	echo "Adjustment factor for MCMC set to ${adjrate}"
	ARGS="${ARGS} -adjrate ${adjrate}"
fi
if [ -z "${d0pi}" ]; then
	echo "Initial delta for the pi all. freq. proposal for MCMC not set"
else 
	echo "Initial delta for the pi all. freq. proposal for MCMC set to ${d0pi}"
	ARGS="${ARGS} -d0pi ${d0pi}"
fi
if [ -z "${upalphaalt}" ]; then
	echo "Alternative update of the pij for MCMC not set"
else 
	echo "Alternative update of the pij for MCMC set"
	ARGS="${ARGS} -upalphaalt ${upalphaalt}"
fi
if [ -z "${d0pij}" ]; then
	echo "Initial delta for the pij all. freq. proposal for MCMC not set"
else 
	echo "Initial delta for the pij all. freq. proposal for MCMC set to ${d0pij}"
	ARGS="${ARGS} -d0pij ${d0pij}"
fi
if [ -z "${d0yij}" ]; then
	echo "Initial delta for the yij all. freq. proposal for MCMC not set"
else 
	echo "Initial delta for the yij all. freq. proposal for MCMC set to ${d0yij}"
	ARGS="${ARGS} -d0yij ${d0yij}"
fi
if [ -z "${seed}" ]; then
	echo "Random Number Generator for MCMC not set"
else 
	echo "Random Number Generator for MCMC set to ${seed}"
	ARGS="${ARGS} -seed ${seed}"
fi




echo "Argument Line:"
echo "i_baypass ${ARGS}"
echo "Starting BayPass"
sources/i_baypass ${ARGS}
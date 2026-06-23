rm(list=ls());

lib_loc="~/R-packages";
library(stargazer,lib.loc=lib_loc);
library(zoo,lib.loc=lib_loc);
library(lmtest,lib.loc=lib_loc);
library(MASS);
library(lme4,lib.loc=lib_loc);
library(arm,lib.loc=lib_loc);
library(base);


### Load data ###
dt=read.table("data/20230723_Mexico_Quarterly.csv",header=T,sep=",");

### load data ###
{
	Forex=dt$Forex;
	GDP=dt$GDP;
	Employment=dt$Employment;
	EmploymentA=dt$EmploymentA;
	EmploymentM=dt$EmploymentM;
	EmploymentS=dt$EmploymentS;
	r_star=dt$r_star;
	r=dt$r;
	bI=dt$bI;
}

### process data ###
{
	forex=Forex/GDP;	# reserve gdp ratio
	### express sectoral employment in per centage of working age population
	Total_Emp=EmploymentA+EmploymentM+EmploymentS;
	empA=(EmploymentA/Total_Emp)*Employment;
	empM=(EmploymentM/Total_Emp)*Employment;
	empS=(EmploymentS/Total_Emp)*Employment;
	emp=Employment/100;
	### express bI as per centage of gdp
	bI=bI/GDP;
}

### adjust data ###
{
	empN=(empA+empS)/100;
	empT=empM/100;

	### r_star and r adjusted to correspond to quarterly values
	r_star=r_star/400;
	r=r/400;
	tau=r-r_star;
}

### normalise data with employment ###
{
	forex=forex*emp;
	bI=bI*emp;
}
### calibrate using normalised data ###
{
	y=forex;
	l=length(y);
	y1=y[1:(l-1)];
	y2=y[2:l];
	lm1=lm(y2~y1);
	names(lm1$coeff)=NULL;
	rho=lm1$coeff[2];
	mu=mean(y);
	sig=sd(y);
	sigma=sig*sqrt(1-rho^2);
	beta=mean(1/(1+r));
	alpha=mean(empN);
	psi=mean(emp);
	Gamma=mean(tau/bI,na.rm=T);

	print(c("beta=",beta),quote=F);
	print(c("alpha=",alpha),quote=F);
	print(c("psi=",psi),quote=F);
	print(c("rho=",rho),quote=F);
	print(c("mu=",mu),quote=F);
	print(c("sigma=",sigma),quote=F);
}

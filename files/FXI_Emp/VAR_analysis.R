rm(list=ls());

lib_loc="~/R-packages";
library(zoo, lib.loc=lib_loc);
library(sandwich, lib.loc=lib_loc);
library(strucchange, lib.loc=lib_loc);
library(urca, lib.loc=lib_loc);
library(lmtest, lib.loc=lib_loc);
library(vars, lib.loc=lib_loc);

dt=read.table("data/20230723_Mexico_Quarterly.csv", header=T, sep=",");

nhor=15;

### selected variables ###
{
	l=dt$Employment/100;
	forex=dt$Forex/dt$GDP;
	gdp=dt$GDP
	df=data.frame(forex,l);
	#df.lev=dt[,c("Forex", "Employment")];
	m=as.matrix(df);
	nr=nrow(df);
}

### VAR model ###
{
	var.model=VAR(df,p=1);
	var.pred=predict(var.model,n.ahead=nhor);
}

### impulse response ###
{
	#feir=irf(var.model, impulse="forex", response="gdp", n.ahead=nhor, ortho=F, runs=1000);
	# plot(feir);
	oir=irf(var.model, impulse="forex", response="l", n.ahead=59, ortho=T, runs=1000, seed=12345);
	# plot(oir);
}

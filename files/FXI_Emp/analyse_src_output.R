rm(list=ls());

source("calibrate.R");
# library(tseries,lib.loc=lib_loc);
library(sandwich,lib.loc=lib_loc);
library(strucchange,lib.loc=lib_loc);
library(urca,lib.loc=lib_loc);
library(lmtest,lib.loc=lib_loc);
library(vars,lib.loc=lib_loc);

### load social planning data ###
{
	Rsp=as.matrix(read.table("output/SP/R.dat", header=F));
	R_starsp=as.matrix(read.table("output/SP/R_star.dat", header=F));
	lNsp=as.matrix(read.table("output/SP/lN.dat", header=F));
	lTsp=as.matrix(read.table("output/SP/lT.dat", header=F));
	lsp=as.matrix(read.table("output/SP/l.dat", header=F));
	bG_starsp=as.matrix(read.table("output/SP/bG_star.dat", header=F));
}

### load competitive equilibrium data ###
{
	Rce=as.matrix(read.table("output/CE/R.dat", header=F));
	R_starce=as.matrix(read.table("output/CE/R_star.dat", header=F));
	lNce=as.matrix(read.table("output/CE/lN.dat", header=F));
	lTce=as.matrix(read.table("output/CE/lT.dat", header=F));
	lce=as.matrix(read.table("output/CE/l.dat", header=F));
}

### data variables already loaded ###
##	r, r_star, empN, empT, emp, forex

### Impulse response function ###
#{
#	forex_sim=(rowSums(bG_starsp,1)/100)[10:59];
#	l_sim=(rowSums(lsp,1)/100)[10:59];
#	df=data.frame(forex_sim,l_sim);
#	m=as.matrix(df);
#	nr=nrow(df);
#	var.model=VAR(df,p=1);
#	oir=irf(var.model, impulse="forex_sim", response="l_sim",n.ahead=40,ortho=T,runs=1000);
#}

### plot and save the data ###
TT=60;
SS=100;

### fig8 ###
{
	png("output/fig/fig8.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(time[2:TT],rowSums(lsp,1)/SS,type='b',pch=15,xlab="Year",ylab="Total employment",ylim=c(0.5,0.75));
	par(new=T);
	lines(time[2:TT],rowSums(lce,1)/SS,type='b',pch=16);
	par(new=T);
	lines(time[2:TT],emp[2:TT],type='b',pch=17);
	legend("topright",legend=c("RP Model","CE Model","Data"),pch=c(15,16,17));
	dev.off();
}
### fig9 ###
{
	png("output/fig/fig9.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(time[2:TT],rowSums(lNsp,1)/SS,type='b',pch=15,xlab="Year",ylab="Non tradable employment",ylim=c(0.3,0.6));
	par(new=T);
	lines(time[2:TT],rowSums(lNce,1)/SS,type='b',pch=16);
	par(new=T);
	lines(time[2:TT],empN[2:TT],type='b',pch=17);
	legend("topright",legend=c("RP Model","CE Model","Data"),pch=c(15,16,17));
	dev.off();
}
### fig10 ###
{
	png("output/fig/fig10.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(time[2:TT],rowSums(lTsp,1)/SS,type='b',pch=15,xlab="Year",ylab="Tradable employment",ylim=c(0,0.4));
	par(new=T);
	lines(time[2:TT],rowSums(lTce,1)/SS,type='b',pch=16);
	par(new=T);
	lines(time[2:TT],empT[2:TT],type='b',pch=17);
	legend("topright",legend=c("RP Model","CE Model","Data"),pch=c(15,16,17));
	dev.off();
}
### fig11 ###
{
	png("output/fig/fig11.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(time[2:TT],rowSums(Rsp,1)/SS,type='b',pch=15,xlab="Year",ylab="Domestic interest rate",ylim=c(1,1.05));
	par(new=T);
	lines(time[2:TT],rowSums(Rce,1)/SS,type='b',pch=16);
	par(new=T);
	lines(time[2:TT],1+r[2:TT],type='b',pch=17);
	legend("topright",legend=c("RP Model","CE Model","Data"),pch=c(15,16,17));
	dev.off();
}
### fig12 ###
{
	png("output/fig/fig12.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(time[2:TT],rowSums(R_starsp,1)/SS,type='b',pch=15,xlab="Year",ylab="World interest rate",ylim=c(1,1.02));
	#par(new=T);
	#lines(time[2:TT],rowSums(R_starce,1)/SS,type='b',pch=16);
	par(new=T);
	lines(time[2:TT],1+r_star[2:TT],type='b',pch=17);
	legend("topright",legend=c("Model","Data"),pch=c(15,17));
	dev.off();
}
### fig13 ###
{
	png("output/fig/fig13.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(time[2:TT],rowSums(bG_starsp,1)/SS,type='b',pch=15,xlab="Year",ylab="Normalised reserve",ylim=c(0,0.5));
	par(new=T);
	lines(time[2:TT],forex[2:TT],type='b',pch=17);
	legend("topright",legend=c("Model","Data"),pch=c(15,17));
	dev.off();
}
### fig14 ###
{
	diff=rowSums(lsp,1)/SS-rowSums(lce,1)/SS;
	diffN=rowSums(lNsp,1)/SS-rowSums(lNce,1)/SS;
	diffL=rowSums(lTsp,1)/SS-rowSums(lTce,1)/SS;
	png("output/fig/fig14.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(time[2:TT],diff,type='b',pch=15,xlab="Year",ylab=expression(Delta));
	legend("topright",legend=c("Emp_RP-Emp_CE"),pch=c(15));
	dev.off();
}

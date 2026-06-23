rm(list=ls());

lib_loc="~/R-packages";
library(stargazer, lib.loc=lib_loc);
library(zoo, lib.loc=lib_loc);
library(lmtest, lib.loc=lib_loc);
library(sandwich, lib.loc=lib_loc);
library(MASS);
library(lme4, lib.loc=lib_loc);
library(arm, lib.loc=lib_loc);

dt=read.table("data/20230121_Mexico_Monthly.csv", header=T, sep=",");

###	read data ###
{
	r=dt$r_star;
	y=dt$forex;
	l=length(y);
}

###	prepare data for AR(1) regression ###
{
	y1=y[1:(l-1)]
	y2=y[2:l];
}

###	perform the regression ###
{
	lm1=lm(y2~y1);
}

###	get the results ###
{
	stargazer(lm1,type="latex");
	res=resid(lm1);	# data on residual is needed to test hogeneity.
	print(bptest(lm1));
}

###	plot and save the data ###
###	fig2 ###
{
	png("output/fig/fig2.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(y1/1000, y2/1000, type='p', pch=15, main="Reserve(bn USD).", xlab=expression(forex[t]), ylab=expression(forex[t+1]), cex.lab=1.5);
	dev.off();
}
###	fig3 ###
{
	png("output/fig/fig3.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(res/1000, type='p', main="Reserve(bn USD).", xlab="Months since Jan 2000", ylab=expression(epsilon[t]), cex.lab=1.5);
	dev.off();
}
###	fig4 ###
{
	png("output/fig/fig4.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	hist(res/1000, breaks=20, freq=T, main="Distribution of residuals (bn USD).", ylab="frequancy", xlab=expression(epsilon[t]), cex.lab=1.5);
	dev.off();
}
###	fig5 ###
{
	png("output/fig/fig5.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	hist(r, breaks=10, freq=T, main="Distribution of world interest rate.", ylab="frequancy", xlab=expression('r*'), cex.lab=1.5);
	dev.off();
}
###	fig6 ###
{
	png("output/fig/fig6.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	discrete.histogram(round(r), freq=T, xlab="r*", prob.col="black")
	dev.off();
}
###	fig7 ###
{
	png("output/fig/fig7.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,5,5,0.2));
	plot(y1/1000,res/1000, type='p', main="", xlab=expression(forex[t]), ylab=expression(epsilon[t]), cex.lab=1.5);
	dev.off();
}

rm(list=ls());

lib_loc="~/R-packages";
library(stargazer,lib.loc=lib_loc);
library(zoo,lib.loc=lib_loc);
library(lmtest,lib.loc=lib_loc);
library(MASS);
library(lme4,lib.loc=lib_loc);
library(arm,lib.loc=lib_loc);
library(base);
library(mFilter,lib.loc=lib_loc);

###	Load data ###
dt1=read.table("data/20230121_Mexico_Monthly.csv",header=T,sep=",");
dt2=read.table("data/20230121_Mexico_Quarterly.csv",header=T,sep=",");
dt3=read.table("data/20230723_Mexico_Quarterly.csv",header=T,sep=",");

### load data ###
{
	Forex=dt3$Forex;
	Population=dt3$Population;
	GDP=dt3$GDP;
	Employment=dt3$Employment;
	EmploymentA=dt3$EmploymentA;
	EmploymentM=dt3$EmploymentM;
	EmploymentS=dt3$EmploymentS;
	r_star=dt3$r_star;
	r=dt3$r;
	bI=dt3$bI;
	year=dt3$Year;
	qtr=dt3$Quarter;
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
	time=as.yearqtr(year+(qtr-1)/4);
}

### adjust data ###
{
	empN=(empA+empS)/100;
	empT=empM/100;
	
	emp.hp=hpfilter(emp, freq=1600, type="lambda");
	empN.hp=hpfilter(empN, freq=1600, type="lambda");
	empT.hp=hpfilter(empT, freq=1600, type="lambda");
	forex.hp=hpfilter(forex, freq=1600, type="lambda");

	emp_cyc=emp.hp$cycle;
	empN_cyc=empN.hp$cycle;
	empT_cyc=empT.hp$cycle;
	forex_cyc=forex.hp$cycle;
	emp_trend=emp.hp$trend;
	empN_trend=empN.hp$trend;
	empT_trend=empT.hp$trend;
	forex_trend=forex.hp$trend;

	### r_star and r adjusted to correspond to quarterly values
	r_star=r_star/400;
	r=r/400;
	tau=r-r_star;
}

### normalise data with employment ###
{
	forex=forex*emp;
	bI=bI*emp;
	forex_cyc=forex_cyc*emp;
	forex_trend=forex_trend*emp;
}
### calibrate rho using nomalised data ###
{
	y=forex;
	l=length(y);
	y1=y[1:(l-1)];
	y2=y[2:l];
	lm1=lm(y2~y1);
	stargazer(lm1,type="text");
	names(lm1$coeff)=NULL;
	rho=lm1$coeff[2];
	mu=mean(forex);
	sig=sd(forex);
	sigma=sig*sqrt(1-rho^2);
	beta=mean(1/(1+r_star));
	Gamma=mean(tau/bI,na.rm=T);
	print(c("beta=",beta));
	print(c("rho=",rho));
	print(c("mu=",mu));
	print(c("sigma=",sigma));
	print(c("Gamma=",Gamma));
}

### caliberate other parameters ###
{
	Gamma=mean(tau/bI,na.rm=T);
	sigma=sd(forex);
	mu=mean(forex);
	beta=mean(1/(1+r_star));
}

### plot data ###
### fig15 ###
{
	png("output/additionalFig/fig15.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,4,5,4));
	plot(time,emp,type='b',pch=15,xlab="Year",ylab="Total Employment");
	par(new=T);
	plot(time,forex,type='b',pch=16, axes=F,xlab="", ylab="");
	axis(side=4,at=pretty(range(forex)));
	mtext("Reserve", side=4, line=3);
	legend("bottomright",legend=c("Employment", "Reserve"), pch=c(15,16));
	dev.off();
}
### fig16 ###
{
	png("output/additionalFig/fig16.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,4,5,4));
	plot(time,empT,type='b',pch=15,xlab="Year",ylab="Tradable Employment");
	par(new=T);
	plot(time,forex,type='b',pch=16,axes=F,xlab="",ylab="");
	axis(side=4,at=pretty(range(forex)));
	mtext("Reserve",side=4, line=3);
	legend("bottomright",legend=c("Tradable Employment","Reserve"), pch=c(15,16));
	dev.off();
}
### fig17 ###
{
	png("output/additionalFig/fig17.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,4,5,4));
	plot(time,empN, type='b', pch=15, xlab="Year", ylab="Non tradable Employment");
	par(new=T);
	plot(time,forex,type='b', pch=16,axes=F,xlab="",ylab="");
	axis(side=4,at=pretty(range(forex)));
	mtext("Reserve", side=4,line=3);
	legend("bottomright", legend=c("Non tradable Employment", "Reserve"), pch=c(15,16));
	dev.off();
}
### fig18 ###
{
	png("output/additionalFig/fig18.png", width=10, height=7.5, units="in", res=200);
	par(mar=c(5,4,5,4));
	plot(time,empN, type='b', pch=15, xlab="Year", ylab="Non tradable Employment");
	par(new=T);
	plot(time,empT,type='b', pch=16,axes=F, xlab="", ylab="");
	axis(side=4,at=pretty(range(empT)));
	mtext("Tradable Employment", side=4, line=3);
	legend("bottomright", legend=c("Non tradable Employment", "Tradable Employment"), pch=c(15,16));
	dev.off();
}

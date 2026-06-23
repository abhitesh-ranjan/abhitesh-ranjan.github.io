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
}

### adjust data ###
{
	forex=Forex/GDP;
	emp=Employment/100;
	forex=forex*emp;
}

### regression ###
{
	y=forex;
	l=length(y);
	y1=y[1:(l-1)];
	y2=y[2:l];
	lm1=lm(y2~y1);
	stargazer(lm1,type="text");
}

rm(list=ls());

lib_loc="~R/packages";
library(stargazer,lib.loc=lib_loc);
library(zoo,lib.loc=lib_loc);
library(lmtest,lib.loc=lib_loc);
library(MASS);
library(lme4,lib.loc=lib_loc);
library(arm,lib.loc=lib_loc);
library(base);
library(mFilter,lib.loc=lib_loc);

### Load data ###
dt=read.table("data/20230121_Mexico_Monthly.csv",header=T,sep=",");

### load data ###
{
	y=dt$forex;
}

l=length(y);
y1=y[1:(l-1)];
y2=y[2:l];
lm1=lm(y2~y1);
stargazer(lm1,type="text");

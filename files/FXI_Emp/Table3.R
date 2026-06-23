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

dt=read.table("data/20230723_Mexico_Quarterly.csv",header=T,sep=",");


### load data ###
{
	Forex=dt$Forex;
	Population=dt$Population;
	GDP=dt$GDP;
	Employment=dt$Employment;
	EmploymentA=dt$EmploymentA;
	EmploymentM=dt$EmploymentM;
	EmploymentS=dt$EmploymentS;
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
}

### adjust data ###
{
	empN=(empA+empS)/100;
	empT=empM/100;
}

### normalise data ###
{
	forex=forex*emp;
}

### display ###
{
	print("======================================",quote=F);
	print("   Variable  | Mean   |   Std.dev.",quote=F);
	print("======================================",quote=F);
	print(c("   Forex     | ",round(mean(forex),3),"|",round(sd(forex),2)),quote=F);
	print(c("  Employment | ",round(mean(emp),3),"|",round(sd(emp),4)),quote=F);
	print(c("Non Tradable Employment |",round(mean(empN),3),"|",round(sd(empN),4)),quote=F);
	print(c("Tradable Employment  | ",round(mean(empT),3),"    | ",round(sd(empT),4)),quote=F);
	print("--------------------------------------",quote=F);
}

rm(list=ls());
lib_loc="~/R-packages";
library(base);
library(plm,lib.loc=lib_loc);
library(stargazer,lib.loc=lib_loc);

dt=read.table("data/data_20230121.csv", header=T, sep=",");
{
	dt$Forex=10^6*dt$Forex;
	dt$FbyG=dt$Forex/dt$GDP;
	#dt$GDP2=dt$GDP^2;
	#dt$dFbyG=0*dt$FbyG;	# initiate dFbyG
	#for (cty in levels(dt$Country))
	#{
	#	for (Y in 2019:2001)
	#	{
	#		dt[dt$Country==cty & dt$Year==Y,]$dFbyG = dt[dt$Country==cty & dt$Year==Y,]$FbyG-dt[dt$Country==cty & dt$Year==(Y-1),]$FbyG
	#	}
	#}
	#dt$CC=1-dt$ka_open;
	#dt1=subset(x=dt[dt$Year!=2000,],subset= (FbyG<1 & !is.na(Forex) & !is.na(GDP) & !is.na(Employment)), select=c(Country,Year,Forex,GDP,Employment, Employment_Ag,FbyG,r_star,CC, GDP2, dFbyG));
	dt1=subset(x=dt, subset= (!is.na(Forex)& !is.na(GDP) & !is.na(GDP) & !is.na(Employment)), select=c(Country,Year,Forex,Employment, FbyG,r_star));
	dt1=pdata.frame(dt1, index=c("Country","Year"))
	#dt1=dt1[complete.cases(dt1),]
	dt1$X=dt1$FbyG*dt1$Employment;
	fit=plm(Employment~X+r_star, data=dt1, model="within", effect="twoways");
	#print(summary(fit));
	#Coeff=summary(fit)$coefficients;
	stargazer(fit,type="text");
}

### calculate FbyG change ###
{
	FbyG2=dt1[dt1$Year==2000,]$FbyG;
	FbyG3=dt1[dt1$Year==2019,]$FbyG;

}
#print(Coeff[2,]);
#for (Y in 2000:2019)
#{
#	df=dt1[dt1$Year==Y,];
#	print(c("Year= ",Y));
#	fit=lm(Employment_Ag~FbyG,data=df);
#	print(summary(fit));
#}
#lag=10;
#for (Y in (2000+lag):2019)
#{
#	if(nrow(dt1[dt1$Year==Y,])>nrow(dt1[dt1$Year==Y-lag,]))
#	{
#		df=merge(dt1[dt1$Year==Y,], dt1[dt1$Year==Y-lag,],by.x="Country",by.y="Country",all.x=F,all.y=T);
#	}
#	else
#	{
#		df=merge(dt1[dt1$Year==Y,], dt1[dt1$Year==Y-lag,],by.x="Country",by.y="Country",all.x=T,all.y=F);
#	}
#	#C=cor.test(df$FbyG,df$Employment);
#	print(c("Year= ",Y));
#	fit=lm(Employment.x~FbyG.x+FbyG.y, data=df, na.action=na.omit);
#	print(summary(fit));
#}
#count=0;
#names=c();
#for (cty in levels(dt1$Country))
#{
#	df=dt1[dt1$Country==cty,];
#	#df=df[complete.cases(df),];
#	print(c("Country= ", cty));
#	if (nrow(df)>1)
#	{
#		fit=lm(Employment~FbyG, data=df, singular.ok=FALSE);
#		Coeff=summary(fit)$coefficients;
#		print(Coeff);
#		if((Coeff[2,1]<0))
#		{
#			count=count+1;
#			names=c(names,cty);
#		}
#	}
#}
#print(c("Count= ",count));
#print(names);

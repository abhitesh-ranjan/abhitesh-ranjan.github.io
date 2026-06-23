rm(list=ls());

lib_loc="~/R-packages";



### load social planning data ###
{
	Vsp=as.matrix(read.table("output/SP/V.dat", header=F));
	Rsp=as.matrix(read.table("output/SP/R.dat", header=F));
	R_starsp=as.matrix(read.table("output/SP/R_star.dat", header=F));
	lNsp=as.matrix(read.table("output/SP/lN.dat", header=F));
	lTsp=as.matrix(read.table("output/SP/lT.dat", header=F));
	lsp=as.matrix(read.table("output/SP/l.dat", header=F));
	bG_starsp=as.matrix(read.table("output/SP/bG_star.dat", header=F));
}

### compute average time series ###
{
	V=rowSums(Vsp)/100;
	R=rowSums(Rsp)/100;
	R_star=rowSums(R_starsp)/100;
	lN=rowSums(lNsp)/100;
	lT=rowSums(lTsp)/100;
	l=rowSums(lsp)/100;
	bG_star=rowSums(bG_starsp)/100;
}

### compute additional variables ###
{
	bG=-bG_star;
	bI=(R-R_star)/1;
	bH=-bG-bI;
}

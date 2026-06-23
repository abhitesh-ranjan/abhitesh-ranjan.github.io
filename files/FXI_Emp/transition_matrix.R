rm(list=ls());

lib_loc="~/R-packages";

dt=read.table("data/20230121_Mexico_Monthly.csv", header=T, sep=",");

###	compute transition matrix ###
{
	r=dt$r_star;
	r=round(r);
	k=length(unique(r));
	M=matrix(0,k,k);
	for(row in 1:k)
	{
		for(col in 1:k)
		{
			count=0;
			i=row-1;
			j=col-1;
			for(m in 2:length(r))
			{
				if((r[m-1]==i)&&(r[m]==j))
				{
					count=count+1;
				}
			}
			M[row,col]=count;
		}
	}
	S=rowSums(M);
	M=M/S;
	M=round(M,4);
	write.matrix(M,file="output/M.dat");
}

/*
In this code, I solve the dynamic stochastic version of two sector economy.
*/

#include <iostream>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <cmath>
#include <string>

using namespace std;
using namespace arma;

arma::mat tauchen(float rho, float mu, float sigma, float m, int n);

class model{
	private:
		// parmaeters
		float beta;	// discount factor
		float alpha;	// trade openess parameter
		float psi;	// labor elasticity
		float Gamma;	// capital control
		float rho;	// persistence parameter
		float mu;	// average reserve level
		float sigma;	// std dev of reserve
		int N;		// number of grid points for domestic rate
		int Sr;		// number of grid points for world rate
		int St;		// number of grid points for government transfer
		mat M;		// transition probability matrix for R*
		mat p;		// transition probability matrix for bG_star
		vec R;
		vec R_star;
		vec bG_star;
		bool SP;	// true for social planner; false for competitive equilibrium

		// solution
		vec w;		// wage rate
		vec R_prime;	// interest rate
		vec bH_prime;	// HH saving
		vec V;		// value function
		vec c;		// consumption
		uvec idx;	// index of optimal R_prime;

		// functions
		bool check_param(void);
		void solve(void);
		mat E(vec X);
		float u(float cc);
		vec u(vec cc);
		mat u(mat cc);
	public:
		model(float bbeta, float aalpha, float ppsi, float GGamma, float rrho, float mmu, float ssigma, int NN, int SSt, mat MM, bool SSP);
		void simulate(int T, int S);
	};

int main(int argc, char** argv)
{
	float beta=0.9;
	float alpha;
	float psi=0.63;
	alpha=1-1.14*psi;
	float Gamma=0.012;
	Gamma=1;
	float rho=0.961;
	float mu=0.295;
	float sigma=0.022;
	int N=20;
	int Sr=8;
	int St=10;
	mat M;
	M.load("output/M.dat",raw_ascii);
	M.print("M=");

	model mdl1(beta,alpha,psi,Gamma,rho,mu,sigma,N,St,M,true);
	mdl1.simulate(59,100);

	return 0;
}

void model::simulate(int T, int S)
{
	float forex0=0.1750412, forex1=0.1749491, r_star0=0.00696, r_star1=0.00808;
	if(!SP)
	{
		forex0=0;
		forex1=0;
	}
	float R_star0=1+r_star0, R_star1=1+r_star1;
	float r0=0.0238;
	float R0=1+r0;
	// Range checks
	bool label=true;
	if((R_star.min()>R_star0)||(R_star.max()<R_star0))
	{
		label=false;
		cerr<<"R_star0 out of range. R_star0="<<R_star0<<" R_star.min="<<R_star.min()<<", R_star.max="<<R_star.max()<<"."<<endl;
	}
	if((R_star.min()>R_star1)||(R_star.max()<R_star1))
	{
		label=false;
		cerr<<"R_star1 out of range. R_star1="<<R_star1<<" R_star.min="<<R_star.min()<<", R_star.max="<<R_star.max()<<"."<<endl;
	}
	if((bG_star.min()>forex0)||(bG_star.max()<forex0))
	{
		label=false;
		cerr<<"forex0 out of range. forex0="<<forex0<<" bG_star.min="<<bG_star.min()<<" bG_star.max="<<bG_star.max()<<"."<<endl;
	}
	if((bG_star.min()>forex1)||(bG_star.max()<forex1))
	{
		label=false;
		cerr<<"forex1 out of range. forex1="<<forex1<<" bG_star.min="<<bG_star.min()<<" bG_star.max="<<bG_star.max()<<"."<<endl;
	}
	if((R.min()>R0)||(R.max()<R0))
	{
		label=false;
		cerr<<"R0 out of range. R0="<<R0<<" R.min="<<R.min()<<" R.max="<<R.max()<<"."<<endl;
	}
	if(!label)
	{
		cout<<"Unable to simulate."<<endl;
		exit(0);
	}

	// Determine initial state
	// R_star
	int sr0, n0, sr1, st0, st1;
	int n1, sr2, st2;
	for(int sr=1;sr<Sr;sr++)
	{
		float R_star_mid=(R_star(sr-1)+R_star(sr))/2;
		if((R_star(sr-1)<=R_star0)&&(R_star0<R_star_mid))
			sr0=sr-1;
		else if((R_star_mid<=R_star0)&&(R_star0<=R_star(sr)))
			sr0=sr;
		if((R_star(sr-1)<=R_star1)&&(R_star1<R_star_mid))
			sr1=sr-1;
		else if((R_star_mid<=R_star1)&&(R_star1<=R_star(sr)))
			sr1=sr;
	}
	// forex
	for(int st=1;st<St;st++)
	{
		float bG_star_mid=(bG_star(st-1)+bG_star(st))/2;
		if((bG_star(st-1)<=forex0)&&(forex0<bG_star_mid))
			st0=st-1;
		else if((bG_star_mid<=forex0)&&(forex0<=bG_star(st)))
			st0=st;
		if((bG_star(st-1)<=forex1)&&(forex1<bG_star_mid))
			st1=st-1;
		else if((bG_star_mid<=forex1)&&(forex1<=bG_star(st)))
			st1=st;
	}
	// R
	for(int n=1;n<N;n++)
	{
		float R_mid=(R(n-1)+R(n))/2;
		if((R(n-1)<=R0)&&(R0<R_mid))
			n0=n-1;
		else if((R_mid<=R0)&&(R0<R(n)))
			n0=n;
	}

	// Declare variables.
	mat sim_idx(T,S);
	int init_state= sr0*N*Sr*St*St+n0*Sr*St*St+sr1*St*St+st0*St+st1;
	int next_state;
	sim_idx.row(0)=init_state*ones<rowvec>(S);
	
	// Simulation codes.
	for(int s=0;s<S;s++)
	{
		init_state=sim_idx(0,s);
		//cout<<"Simulating column: "<<s<<endl;
		bool valid=false;
		int iter=0;
		do{
			iter++;
			if(iter>1)
			{
				cout<<"Simulating column: "<<s;
				cout<<" Iteration: "<<iter<<endl;
			}
			vec Ur=randu<vec>(T);
			vec Ut=randu<vec>(T);
			for(int t=1;t<T;t++)
			{
				float sumr0=0, sumr1=0;
				// Get values for sr2, st2
				for(int j=0;j<Sr;j++)
				{
					sumr1=sumr0+M(sr1,j);
					if((sumr0<Ur(t))&&(Ur(t)<=sumr1))
						sr2=j;
					sumr0=sumr1;
				}
				float sumt0=0, sumt1=0;
				for(int j=0;j<St;j++)
				{
					sumt1=sumt0+p(st1,j);
					if((sumt0<Ut(t))&&(Ut(t)<=sumt1))
						st2=j;
					sumt0=sumt1;
				}
				n1=idx(init_state);
				next_state=sr1*N*Sr*St*St+n1*Sr*St*St+sr2*St*St+st1*St+st2;
				sim_idx(t,s)=next_state;
				init_state=next_state;
			}
			// Check validity
			vec Vcheck(T);
			for(int t=0;t<T;t++)
				Vcheck(t)=V(sim_idx(t,s));
			if(Vcheck.is_finite())
				valid=true;
			//valid=true;
		}while(!valid);
	}
	//cout<<sim_idx;
	mat Vmat(T,S);
	mat Rmat(T,S);
	mat R_mat(T,S);
	mat lNmat(T,S), lTmat(T,S), lmat(T,S);
	mat forex_mat(T,S);
	int count=0;
	for(int t=0;t<T;t++)
	{
		for(int s=0;s<S;s++)
		{
			Vmat(t,s)=V(sim_idx(t,s));
			Rmat(t,s)=R_prime(sim_idx(t,s));
			int r=sim_idx(t,s);
			int sr1=(r/(St*St))%Sr;
			R_mat(t,s)=R_star(sr1);
			lmat(t,s)=1-psi*c(sim_idx(t,s));
			lNmat(t,s)=(1-alpha)*c(sim_idx(t,s));
			if(lmat(t,s)<=0)
			{
				count++;
				cerr<<count<<": Negative employment."<<endl;
			}
			int st=(r/St)%St;
			int st1=r%St;
			forex_mat(t,s)=bG_star(st1);
		}
	}
	lTmat=lmat-lNmat;
	uvec non_fin=find_nonfinite(Vmat);
	cout<<"# non finite values in Vmat: "<<non_fin.n_elem<<endl;
	//cout<<"Rmat: "<<Rmat;

	// Compute required moments
	float lmean=mean(mean(lmat));
	float lvar=mean(var(lmat));
	float lNmean=mean(mean(lNmat));
	float lNvar=mean(var(lNmat));
	float lTmean=mean(mean(lTmat));
	float lTvar=mean(var(lTmat));
	cout<<"Employment"<<endl;
	cout<<"Mean: "<<lmean<<" Var: "<<lvar<<endl;
	cout<<"Non tradable employment"<<endl;
	cout<<"Mean: "<<lNmean<<" Var: "<<lNvar<<endl;
	cout<<"Tradable employment"<<endl;
	cout<<"Mean: "<<lTmean<<" Var: "<<lTvar<<endl;
	cout<<"Domestic rate"<<endl;
	cout<<"Mean: "<<mean(mean(Rmat))<<endl;
	cout<<"International rate"<<endl;
	cout<<"Mean: "<<mean(mean(R_mat))<<endl;

	// Save data
	string loc;
	if(SP)
		loc="output/SP/";
	else
		loc="output/CE/";
	Vmat.save(loc+"V.dat",raw_ascii);
	Rmat.save(loc+"R.dat",raw_ascii);
	R_mat.save(loc+"R_star.dat",raw_ascii);
	lNmat.save(loc+"lN.dat",raw_ascii);
	lTmat.save(loc+"lT.dat",raw_ascii);
	lmat.save(loc+"l.dat",raw_ascii);
	forex_mat.save(loc+"bG_star.dat",raw_ascii);
}
arma::mat model::E(arma::vec X)
{
	int ROW=Sr*N*Sr*St*St;
	int COL=N;
	mat EX=zeros<mat>(ROW,COL);
	for(int r=0;r<ROW;r++)
	{
		int sr=r/(N*Sr*St*St);
		int n=(r/(Sr*St*St))%N;
		int sr1=(r/(St*St))%Sr;
		int st=(r/St)%St;
		int st1=r%St;
		for(int c=0;c<COL;c++)
		{
			int n1=c;
			for(int sr2=0;sr2<Sr;sr2++)
				for(int st2=0;st2<St;st2++)
				{
					int r1=sr1*N*Sr*St*St+n1*Sr*St*St+sr2*St*St+st1*St+st2;
					EX(r,c)=EX(r,c)+p(st1,st2)*M(sr1,sr2)*X(r1);
				}
		}
	}
	return EX;
}

void model::solve(void)
{
	int ROW=Sr*N*Sr*St*St;
	int COL=N;
	R=linspace(1.00,1.05,N);
	R_star=linspace(1.00,1.0175,Sr);
	float sig=sigma/pow(1-rho*rho,0.5);
	int m=2;
	bG_star=linspace(mu-m*sig,mu+m*sig,St);
	if(!SP)
		bG_star=0*bG_star;
	R_prime=R(0)*ones<vec>(ROW);

	V=zeros<vec>(ROW);
	vec V_next=zeros<vec>(ROW);
	mat V_grid(ROW,COL);

	c=zeros<vec>(ROW);
	mat c_grid(ROW,COL);
	bH_prime=zeros<vec>(ROW);

	idx=zeros<uvec>(ROW);
	float err=10, tol=1e-3;
	int iter=0, iter_max=1e6;

	do{
		for(int row=0;row<ROW;row++)
		{
			int sr=row/(N*Sr*St*St);
			int n=(row/(Sr*St*St))%N;
			int sr1=(row/(St*St))%Sr;
			int st=(row/St)%St;
			int st1=row%St;
			for(int col=0;col<COL;col++)
			{
				c_grid(row,col)=1+R_star(sr)*(bG_star(st)-1/Gamma*(R(n)-R_star(sr)))-(bG_star(st1)-1/Gamma*(R(col)-R_star(sr1)))-1/Gamma*pow(R(n)-R_star(sr),2);
			}
		}
		c_grid=c_grid/(1+psi);
		V_grid=u(c_grid)+beta*E(V);
		idx=index_max(V_grid,1);
		for(int r=0;r<ROW;r++)
		{
			V_next(r)=V_grid(r,idx(r));
			R_prime(r)=R(idx(r));
			c(r)=c_grid(r,idx(r));
		}
		//cout<<datum::inf-(datum::inf)<<endl;
		err=abs(V_next-V).max();
		V=V_next;
		//cout<<find_nonfinite(V)<<endl;
		if(iter%100==0)
			cout<<iter<<": "<<err<<endl;
		//cout<<"V="<<V<<endl;
		iter++;
	}while((err>tol)&&(iter<iter_max));
	uvec non_fin=find_nonfinite(V);
	uvec unik=find_unique(idx);
	// count valid 'c'
	int count=sum(psi*c>=1);	
	//cout<<non_fin<<endl;
	//cout<<unik;
	cout<<"# invalid c: "<<count<<endl;
	cout<<"# unique vals: "<<unik.n_elem<<endl;
	cout<<"# non finite: "<<non_fin.n_elem<<endl;
	cout<<iter<<": "<<err<<endl;
}

float model::u(float cc)
{
	if((cc>0))
		return log(cc);
	else
		return -datum::inf;
}
arma::vec model::u(arma::vec cc)
{
	int L=cc.n_elem;
	arma::vec uu(L);
	for(int l=0;l<L;l++)
		uu(l)=u((float)cc(l));
	return uu;
}
arma::mat model::u(arma::mat cc)
{
	int col=cc.n_cols;
	int row=cc.n_rows;
	arma::mat uu(row,col);
	for(int c=0;c<col;c++)
		uu.col(c)=u((arma::vec)cc.col(c));
	return uu;
}

model::model(float bbeta, float aalpha, float ppsi, float GGamma, float rrho, float mmu, float ssigma, int NN, int SSt, mat MM, bool SSP)
{
	SP=SSP;
	if(SP)
		cout<<"___________ Ramsey Planner Problem _______"<<endl;
	else
		cout<<"___________ Competitive Equilibrium _______"<<endl;
	beta=bbeta;
	alpha=aalpha;
	psi=ppsi;
	Gamma=GGamma;
	rho=rrho;
	mu=mmu;
	sigma=ssigma;
	N=NN;
	//Sr=SSr;
	St=SSt;
	M=MM;
	Sr=M.n_rows;
	int m=2;
	p=tauchen(rho,mu,sigma,m,St);
	if(!check_param())
		exit(0);
	solve();
}
bool model::check_param(void)
{
	bool label=true;
	if((beta<=0)||(beta>=1))
	{
		label=false;
		cerr<<"Invalid beta: "<<beta<<endl;
	}
	if((alpha<=0)||(alpha>=1))
	{
		label=false;
		cerr<<"Invalid alpha: "<<alpha<<endl;
	}
	if(psi<0)
	{
		label=false;
		cerr<<"Invalid psi: "<<psi<<endl;
	}
	if(Gamma<=0)
	{
		label=false;
		cerr<<"Invalid Gamma: "<<Gamma<<endl;
	}
	if(abs(rho)>=1)
	{
		label=false;
		cerr<<"Invalid rho: "<<rho<<endl;
	}
	if(mu<=0)
	{
		label=false;
		cerr<<"Invalid mu: "<<mu<<endl;
	}
	if(sigma<0)
	{
		label=false;
		cerr<<"Invalid sigma: "<<sigma<<endl;
	}
	if((N<2)||(Sr<2)||(St<2))
	{
		label=false;
		cerr<<"Invalid grid size. N= "<<N<<", Sr= "<<Sr<<", St= "<<St<<endl;
	}
	if((M.n_rows!=M.n_cols)||any(vectorise(M)<0)||any(vectorise(M)>1)||any(abs(sum(M,1)-1)>1e-3))
	{
		label=false;
		cerr<<"Invalid transition matrix. M="<<M<<endl;
	}
	return label;
}
arma::mat tauchen(float rho, float mu, float sigma, float m, int n)
{
	bool label=true;
	string str="In arma::mat tauchen(float rho, float mu, float sigma, float m, int n)";
	if(abs(rho)>=1)
	{
		label=false;
		cerr<<str<<", invalid rho: "<<rho<<endl;
	}
	if(sigma<0)
	{
		label=false;
		cerr<<str<<", invalid sigma: "<<sigma<<endl;
	}
	if(m<=0)
	{
		label=false;
		cerr<<str<<", invalid m: "<<m<<endl;
	}
	if(n<2)
	{
		label=false;
		cerr<<str<<", invalid n: "<<n<<endl;
	}
	if(!label)
		exit(0);
	float sig=sigma/pow(1-rho*rho,0.5);
	vec x=linspace(mu-m*sig, mu+m*sig,n);
	float w=x(1)-x(0);
	mat p(n,n);
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
		{
			if(j==0)
				p(i,j)=normcdf((x(j)+w/2-(1-rho)*mu-rho*x(i))/sig);
			else if(j==(n-1))
				p(i,j)=1-normcdf((x(j)-w/2-(1-rho)*mu-rho*x(i))/sig);
			else
				p(i,j)=normcdf((x(j)+w/2-(1-rho)*mu-rho*x(i))/sigma)-normcdf((x(j)-w/2-(1-rho)*mu-rho*x(i))/sigma);
		}
	return p;
}

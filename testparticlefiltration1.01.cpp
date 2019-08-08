#include "./lbm/Domain.h"

struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    double rhos;
    double rho1;
    double rho2;
//    double pl;
    double sy;
    int pnx;
    int pny;
    double ratiol;
    double ratiot;
	double ratiom;
	size_t ny;
	double l;
	
};

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}

void Setup(LBM::Domain &dom, void *UD)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    myUserData &dat = (*static_cast<myUserData *> (UD));

    // std::cout<<vvb<<std::endl;
    //specific top layer pressure
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0;ix<nx;++ix)
    {
        double *f = dom.F[ix][ny-1][0];
        double *f1 = dom.F[ix][ny-2][0];
        double rho1 = dom.Rho[ix][ny-2][0];
        Vec3_t vel1 = dom.Vel[ix][ny-2][0];
        
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,dat.rho1,vel1) + f1[k] - dom.Feq(k,rho1,vel1);
        }
        Vec3_t idx(ix,ny-1,0);
        dom.CalcPropsForCell(idx);
    }

    //specific bottom layer pressure
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0;ix<nx;++ix)
    {
        double *f = dom.F[ix][0][0];
        double *f1 = dom.F[ix][1][0];
        double rho1 = dom.Rho[ix][1][0];
        Vec3_t vel1 = dom.Vel[ix][1][0];
        
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,dat.rho2,vel1) + f1[k] - dom.Feq(k,rho1,vel1);
        }
        Vec3_t idx(ix,0,0);
        dom.CalcPropsForCell(idx);
    }

    //fix the moving particle if leave the domain in y direction
	int ipp = -1;
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        if(!dom.Particles[ip].IsFree()) continue;
        if(dom.Particles[ip].X(1)<2*(double) dat.R)
        {
            dom.Particles[ip].V = 0.0;
            dom.Particles[ip].W = 0.0;
			dom.Particles[ip].X(1) = -100;
			dom.Particles[ip].Xb(1) = -100;
            dom.Particles[ip].FixVeloc();
        }
		if(dom.Particles[ip].X(1)>ny-1-3*dat.R-dat.l)
        {
            #pragma omp atomic
            ipp++;
		}
	}

    // std::cout<<ipp<<std::endl;
	if(ipp<0)
    {
        double py = dat.ny -1-dat.R;
        Vec3_t pos(0,0,0);
        Vec3_t v(0,0,0);
        Vec3_t w(0,0,0);
        int pnum = dom.Particles.size();
        for(int i=0; i<dat.pnx; ++i)
        {
            // Vec3_t dxr(0.0,0.0,0.0);
			Vec3_t dxr(random(-0.5*dat.R,0.5*dat.R),random(-0*dat.R,0*dat.R),0.0);
            double px = 0.5*dat.l+2*i*dat.R+i*dat.l;
            pos = px, py , 0;
            dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, dat.rhos, dat.R, dom.dtdem));
            // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
            pnum++;
            dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
            // dom.Particles.back().Ff = 0.0, 0.0, 0.0;
            dom.Particles.back().Kn = 50;
            dom.Particles.back().Gn = -0.3;
            dom.Particles.back().Kt = 0.0;
            dom.Particles.back().Mu = 0.4;
            dom.Particles.back().Eta = 0.0;
            dom.Particles.back().Beta = 0.0;
            dom.Particles.back().A = 3e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
            dom.Particles.back().kappa = 1e9*dat.ratiol;
            dom.Particles.back().Z = 3.97e-11*dat.ratiot*dat.ratiot/(dat.ratiol*dat.ratiom);
            dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
            dom.Particles.back().D = 2;
            dom.Particles.back().Rh = 0.80*dat.R;
        }
    
    }

}
void Initial(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        Vec3_t vtemp(0, 0, 0);
        dom.Rho[ix][iy][0] = 1.0;
        dom.Vel[ix][iy][0] = vtemp;
        dom.BForce[ix][iy][0] = 0.0, 0.0, 0.0;
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            dom.F[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
            dom.Ftemp[ix][iy][0][k] = dom.Feq(k,1.0,vtemp);            
        }
    }
}
///////////////////////////////////////////////////////////////////
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL)); 
    size_t Nproc = 5;
    if(argc>=2) Nproc = atoi(argv[1]);
    //slurry///////////////////////////////////////////////
    double nu   = 0.10;
    double R    = 10.0;
    int    pnx  = 10;                    //small particle number
    int    pny  = 5;
	double l    = 19.831;                //小颗粒之间表面间距
    double rho  = 1.0;                   //所有计算运用的是密度比，所以密度的单位不用转换。
    double rhos = 2.7;
    double Ga   = 0;                     // Immersed boundary lattice boltzmann simulation of turbulent channel flows in the presence of spherical particles. Numerical Simulation of Two Spheres with Different Density Settling In Fluids through Lattice Boltzmann Method, 文章中定义了伽利略数的计算公式
    double gy   = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    //filter///////////////////////////////////////////////
    double RR   = 20.0;
    int    Pnx  = 4;                     //big particle number
    int    Pny  = 3;             
    double pdx  = 60;                    //gap between big particle
    double pdy  = 60;
    ////////////////////////////////////////////////////// 
    size_t nx = std::max(std::ceil(Pnx*(RR*2+pdx)), std::ceil(pnx*(R*2+l)));	
	size_t ny = std::ceil((Pny-1)*(std::sqrt(3)*RR+pdy) + 40.0+2*RR + pny*2*R+pny*l)+1;
    size_t nz = 1;	
    double dx = 1.0;
    double dt = 1.0;
    double ratiol = 1e-6; //Lr/Ll
    double ratiot = 1e-7; //Tr/Tl 
	double ratiom = 1e-15; //Mr/Ml    	
	double sy = (Pny-1)*(std::sqrt(3)*RR+pdy) + 2*RR + 0.5*l + R + 40 + 1;
	////////output//////////////////////////////////////////////////////
	std::cout<<"1 water content = "<<(((l+2*R)*(l+2*R)-3.14*R*R)/(3.14*R*R*2.7))<<std::endl;
	std::cout<<"2 AOS/D = "<<pdx/(2*R)<<std::endl;
	std::cout<<"3 filter = "<<Pny<<std::endl;
	std::cout<<"4 R = "<<R<<std::endl;
    std::cout<<"5 RR = "<<RR<<std::endl;
	std::cout<<"nx1="<<std::ceil(Pnx*(RR*2+pdx))<<"; nx2="<<std::ceil(pnx*(R*2+l))<<std::endl;
	std::cout<<"6 nx="<<nx<<"; ny="<<ny<<"; nz="<<nz<<std::endl;
	std::cout<<"7 gy = "<<gy<<std::endl;
	//////////////////////////////////////////////////////////////////////
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);         //去Domain.h文件看一段语句。
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = gy;
    my_dat.R = R;
    my_dat.rho1 = 1.03;
    my_dat.rho2 = 0.97;
	std::cout<<"8 pressure drop = "<<100*(my_dat.rho1-my_dat.rho2)/3<<"kPa"<<std::endl;
    my_dat.pnx = pnx;
    my_dat.pny = pny;
    my_dat.sy = sy;
	my_dat.ny = ny;
	my_dat.l = l;
    my_dat.rhos = rhos;
	my_dat.ratiol = ratiol;
	my_dat.ratiom = ratiom;
	my_dat.ratiot = ratiot;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;  
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;
	//////////////////////////////////////////////////////////////////////
	//////breakpoint executing///////////////////////////////////////////
	Initial(dom,dom.UserData);
    bool iscontinue = false;
	// dom.InitialFromH5("test_wu2_0072.h5",g0);
    // dom.Initial(rho,v0,g0);
    // bool iscontinue = true;
	/////////////draw the filter and particles in slurry//////////////////////////////////////////
    if(!iscontinue)
    {	
		//filter position//////////////////////////////////////////////////////////////////////////
		Vec3_t pos(0.0,0.0,0.0);
		Vec3_t pos1(0.0,0.0,0.0);
		Vec3_t dxp(0.0,0.0,0.0);
		Vec3_t v(0.0,0.0,0.0);
		Vec3_t w(0.0,0.0,0.0);
		double py = 0;
		double px = 0;
		int pnum = 0;
		for(int j=0; j<Pny; ++j)
		{
			py = 40.0 + j*(std::sqrt(3)*RR+pdy) + RR + 1;
			int temp = j%2==0 ? Pnx : Pnx+1;
			for(int i = 0; i<temp; ++i)
			{
				px = j%2==0 ? (0.5*pdx+RR+i*(2*RR+pdx)) : i*(2*RR+pdx);
				pos = px, py, 0;
				dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, RR, dom.dtdem));
				dom.Particles[pnum].FixVeloc();
				pnum++;
			}
		}
		//particle in slurry position/////////////////////////////////////////////////////////////////
		for(int j=0; j<pny; ++j)
		{
			py = sy + j*(2*R+l);
			for(int i=0; i<pnx; ++i)
			{
				Vec3_t dxr(random(-0.5*R,0.5*R),random(-0.5*R,0.5*R),0.0);
				// Vec3_t dxr(0.0,0.0,0.0);
				px = 0.5*l+i*(2*R+l);
				pos = px, py , 0;
				dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
				// std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
				pnum++;
			}
		}
    }
    std::cout<<"9 Particles number = "<<dom.Particles.size()<<std::endl;
    ///////define the properties of filter and particles in the slurry///////////////////////////////////////
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Kn        =  50;
        dom.Particles[ip].Gn        = -0.3;
        dom.Particles[ip].Kt        =  0.0;
        dom.Particles[ip].Mu        =  0.4;
        dom.Particles[ip].Eta       =  0.0;
        dom.Particles[ip].Beta      =  0.0;
        dom.Particles[ip].A         =  3e-20/(((ratiol/ratiot)*(ratiol/ratiot))*ratiom);
        dom.Particles[ip].kappa     =  1e9*ratiol;
        dom.Particles[ip].Z         =  3.97e-11*ratiot*ratiot/(ratiol*ratiom);
        dom.Particles[ip].VdwCutoff =  std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        dom.Particles[ip].D         =  2;
		dom.Particles[ip].Ff        =  0.0, -M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*my_dat.g,0.0, 0.0;
        if(dom.Particles[ip].IsFree())
        {  
            dom.Particles[ip].Rh = 0.80*R;
        }else
        {
            dom.Particles[ip].Rh = 0.80*RR;
        }
    }
	//std::cout<<"VdwCutoff "<<dom.Particles.back().VdwCutoff<<std::endl;
    //OUTPUT set////////////////////////////////////////////////////////////////////////////////////////////////
    double Tf = 1000;
	//std::cout<<"9 Tf="<<Tf<<std::endl;
    dom.IsF = true;    
    double dtout = 100;
    // periodic in x
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    // periodic in y
    // dom.Box = 0.0, ny-1, 0.0;
    // dom.modexy = 1;
    //solving
    dom.SolveIBM( Tf, dtout, "testAN5", Setup, NULL);
}MECHSYS_CATCH
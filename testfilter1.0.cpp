//##################################################################################################################################################################################################
//This is for DATA analysis in the filtration simulation test.
//Special for the initialize the filter.
//For testing testfilter.cpp.
//It is strange that the permeability and hydraulitic coeffient of the filter will change with the pressure drop, but not too much. 
//So in testing the permeability of filter, we fix the pressure drop with 0.5 kPa.
//Built by Silin WU from Hohai University, Nanjing, China. Email: wusilinhhu@126.com
//Before test you only need to modify line 111-114
//The test report can be seen in a word named "10.Filter shentouxishu de wucha ceshi baogao" in chinese.
//BE HAPPY WITH THIS DATA ANALYSIS.   -----data update:2019-07-18
//##################################################################################################################################################################################################

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
    double ratiol;
    double ratiot;
	double ratiom;
	size_t ny;
	
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

int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL)); 
	printf("\n%s--- test for caculate fliter permeabilty 基础条件 ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);	
    size_t Nproc = 5;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // the parameters effect the filter permeability are RR, Pnx, Pny, pdx, pdy , in my opinion, you only need to change these four parameters in each test////////////////////////////////////////////////////////////////////////////////
    double RR = 20;
    int Pny = 3; 
    double pdx = 117;//gap between big particle
    double pdy = 117;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////// don't effect fiter permeability///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double R = 10;
    double nu = 0.10;
    int Pnx = 1;//big particle number
    int pny = 3;
    double l = 3;
    if(argc>=2) Nproc = atoi(argv[1]); 
    size_t nx = std::ceil(Pnx*(RR*2+pdx));	
	size_t ny = std::ceil((Pny-1)*(std::sqrt(3)*RR+pdy) + 20.0+2*RR + pny*2*R+pny*l)+1;
    size_t nz = 1;	
    double dx = 1.0;
    double dt = 1.0;
    double ratiol = 1e-6; //Lr/Ll
    double ratiot = 1e-7; //Tr/Tl 
	double ratiom = 1e-15; //Mr/Ml    	
    double rho = 1.0; 
    double rhos = 2.7;
    double Ga = 0.0; // Immersed boundary lattice boltzmann simulation of turbulent channel flows in the presence of spherical particles. Numerical Simulation of Two Spheres with Different Density Settling In Fluids through Lattice Boltzmann Method, 文章中定义了伽利略数的计算公式
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
	double sy = (Pny-1)*(std::sqrt(3)*RR+pdy) + 2*RR + 0.5*l + R + 20;
	
    std::cout<<"Pnx"<<Pnx<<std::endl;
	std::cout<<"Pny"<<Pny<<std::endl;
    std::cout<<"pdx"<<pdx<<std::endl;
    std::cout<<"pdy"<<pdy<<std::endl;
    std::cout<<"RR = "<<RR<<std::endl;
	std::cout<<"nx="<<std::ceil(Pnx*(RR*2+pdx))<<std::endl;
	std::cout<<"ny="<<std::ceil((Pny-1)*(std::sqrt(3)*RR+pdy) + 20.0+2*RR + pny*2*R+pny*l)+1<<std::endl;
	
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);         //去Domain.h文件看一段语句。
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = gy;
    my_dat.R = R;
    my_dat.rho1 = 1.0075;
    my_dat.rho2 = 0.9925;
	std::cout<<"pressure drop="<<100*(my_dat.rho1-my_dat.rho2)/3<<"kPa"<<std::endl;
    my_dat.sy = sy;
	my_dat.ny = ny;
    my_dat.rhos = rhos;
	my_dat.ratiol = ratiol;
	my_dat.ratiom = ratiom;
	my_dat.ratiot = ratiot;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;  
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;
	
	
	Initial(dom,dom.UserData);
    bool iscontinue = false;
	// dom.InitialFromH5("test_wu2_0072.h5",g0);
    // dom.Initial(rho,v0,g0);
    //bool iscontinue = true;
	
    if(!iscontinue)
    {	
		//fix
		Vec3_t pos(0.0,0.0,0.0);
		Vec3_t v(0.0,0.0,0.0);
		Vec3_t w(0.0,0.0,0.0);
		double py = 0;
		double px = 0;
		int pnum = 0;
		for(int j=0; j<Pny; ++j)
		{
			py = 20.0 + j*(std::sqrt(3)*RR+pdy) + RR + 1;
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
    }

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;

    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Kn = 55.0;
        dom.Particles[ip].Gn = 0.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].A = 3e-20/(((ratiol/ratiot)*(ratiol/ratiot))*ratiom);
        dom.Particles[ip].kappa = 1e9*ratiol;
        dom.Particles[ip].Z = 3.97e-11*ratiot*ratiot/(ratiol*ratiom);
//		dom.Particles[ip].bbeta = 0.3;
//		dom.Particles[ip].epsilon = 2.05769e-20/((ratiol/ratiot)*(ratiol/ratiot));
//		dom.Particles[ip].s = 200e-9/ratiol;
//		dom.Particles[ip].Lc = 100e-9/ratiol;
//		dom.Particles[ip].l = 3.04e-10/ratiol;
        dom.Particles[ip].VdwCutoff = std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        dom.Particles[ip].D = 2;
		dom.Particles[ip].Ff =  0.0, -M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*my_dat.g,0.0, 0.0;
        if(dom.Particles[ip].IsFree())
        {
            
            dom.Particles[ip].Rh = 0.80*R;

        }else{
            dom.Particles[ip].Rh = 0.80*RR;

        }
        


    }
    
	std::cout<<"VdwCutoff "<<dom.Particles.back().VdwCutoff<<std::endl;


    double Tf = 1e5;
	std::cout<<"Tf="<<Tf<<std::endl;
    dom.IsF = true;    
    double dtout = 1e2;
    // periodic in x
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    // periodic in y
    // dom.Box = 0.0, ny-1, 0.0;
    // dom.modexy = 1;
    
    //solving
    dom.SolveIBM( Tf, dtout, "testfilter", Setup, NULL);
    
}MECHSYS_CATCH

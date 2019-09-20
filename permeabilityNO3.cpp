#include "./lbm/Domain.h"

struct myUserData
{
    std::ofstream oss_ss;
    double nu;
    double rhos;
    double rho1;
    double rho2;
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
///////////////////////////////////////////////////////////////////
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL)); 
    size_t Nproc = 5;
    if(argc>=2) Nproc = atoi(argv[1]);
    //slurry///////////////////////////////////////////////
    double nu   = 0.10;
    double rho  = 1.0;                   //所有计算运用的是密度比，所以密度的单位不用转换。
    double rhos = 2.7;
    //double Ga   = 20;                     // Immersed boundary lattice boltzmann simulation of turbulent channel flows in the presence of spherical particles. Numerical Simulation of Two Spheres with Different Density Settling In Fluids through Lattice Boltzmann Method, 文章中定义了伽利略数的计算公式
    ////////////////////////////////////////////////////// 
    size_t nx = 500;	
	size_t ny = 350;
    size_t nz = 1;	
    double dx = 1.0;
    double dt = 1.0;
    double ratiol = 1e-6; //Lr/Ll
    double ratiot = 1e-7; //Tr/Tl 
	double ratiom = 1e-15; //Mr/Ml    	
	//////////////////////////////////////////////////////////////////////
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);         //去Domain.h文件看一段语句。
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.rho1 = 1.03;                                               //******************
    my_dat.rho2 = 0.97;                                               //******************
	std::cout<<"8 pressure drop = "<<100*(my_dat.rho1-my_dat.rho2)/3<<"kPa"<<std::endl;
	my_dat.ny = ny;
    my_dat.rhos = rhos;
	my_dat.ratiol = ratiol;
	my_dat.ratiom = ratiom;
	my_dat.ratiot = ratiot;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;  
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
        //read
        std::fstream ifile("NO.3.txt",std::ios::in);
        int N = 0;
        int pnum = 0;
        Vec3_t v(0.0,0.0,0.0);
		Vec3_t w(0.0,0.0,0.0);
        Vec3_t pos(0.0,0.0,0.0);
        if(!ifile.fail())
        {
            ifile>>N;
            std::cout<<N<<std::endl;
            for(size_t i=0; i<N; ++i)
            {
                double x,y,r;
                ifile>>x>>y>>r;
                pos = x, y + 50, 0;
                dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, r, dom.dtdem));
				dom.Particles[pnum].FixVeloc();
                dom.Particles[pnum].Rh = 0.80*r;
                pnum++;
            }
        }
    }
    //OUTPUT set////////////////////////////////////////////////////////////////////////////////////////////////
    double Tf = 5;                                               //******************
	//std::cout<<"9 Tf="<<Tf<<std::endl;
    dom.IsF = true;    
    double dtout = 1;                                               //******************
    // periodic in x
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    // periodic in y
    // dom.Box = 0.0, ny-1, 0.0;
    // dom.modexy = 1;
    //solving
    std::cout<<1<<std::endl;

    dom.SolveIBM( Tf, dtout, "NO.3", Setup, NULL);                                                 //******************
}MECHSYS_CATCH
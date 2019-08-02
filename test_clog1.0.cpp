#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    double rhos;
    double vb;
};
void Setup(LBM::Domain &dom, void *UD)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    myUserData &dat = (*static_cast<myUserData *> (UD));
    Vec3_t vvb(dat.vb,0.0,0.0);
    // std::cout<<vvb<<std::endl;
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0;ix<nx;++ix)
    {
        double *f = dom.F[ix][ny-1][0];
        double *f1 = dom.F[ix][ny-2][0];
        double rho1 = dom.Rho[ix][ny-2][0];
        Vec3_t vel1 = dom.Vel[ix][ny-2][0];
        
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,rho1,vvb) + f1[k] - dom.Feq(k,rho1,vel1);
        }
        Vec3_t idx(ix,ny-1,0);
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
        Vec3_t vtemp((double) dat.vb*iy/(ny-1), 0, 0);
        // Vec3_t vtemp((double) dat.vb, 0, 0);
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

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 5;
    if(argc>=2) Nproc = atoi(argv[1]); 

    double nu = 0.1;
    double R = 10.0;
    double RR = 10000.0;
    double vb = 0.4; 
    

    //size_t nx = std::ceil(Pnx*(RR*2 + pdx));
    size_t nx = 50;
    //size_t ny = std::ceil((Pny-1)*(std::sqrt(3)*RR+pdy) + 2*RR + ppl + (pny-1)*pl)+6*Rn + 1;
    size_t ny = 60;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    // reality convert to lattice
    double ratiol = 1e-6;                       //Lr/Ll
    double ratiot = 1e-7;                       //Tr/Tl 
	double ratiom = 1e-15;                      //Mr/Ml  
    //std::cout<<"R = "<<R<<std::endl;
    //std::cout<<"RR = "<<RR<<std::endl;
    double rho = 1.0;
    double rhos = 2.7;
    double Ga = 0.0;
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    //std::cout<<"gy = "<<gy<<std::endl;
    //nu = 1.0/30.0;
    //std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = gy;
    my_dat.R = R;
    my_dat.vb = vb;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;       
    
    for(size_t ix=0; ix<nx; ++ix)
    {
        dom.IsSolid[ix][0][0] = true;
    }
    //dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
   
    //initial
    // dom.InitialFromH5("test_cong_0999.h5",g0);
    // dom.Initial(rho,v0,g0);
    Initial(dom,dom.UserData);
    bool iscontinue = false;
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;

    if(!iscontinue)
    {
        //fix
        Vec3_t pos(0.0,0.0,0.0);
        Vec3_t v(0.0,0.0,0.0);
        Vec3_t w(0.0,0.0,0.0);
        int pnum = 0;
        double py = -9980;   
        double px = 25;
        pos = px, py, 0;
        dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, RR, dom.dtdem));
        dom.Particles[pnum].FixVeloc();
        pnum++;    
    
        //move
        py = 30;
        px = 25;
        pos = px, py , 0;
        dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
        // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
        pnum++;
    }
    
    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Kn = 50.0;
        dom.Particles[ip].Gn = 0.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].A = 2.35e-21/(((ratiol/ratiot)*(ratiol/ratiot))*ratiom);
        dom.Particles[ip].kappa = 1e9*ratiol;
        dom.Particles[ip].Z = 3.97e-11*ratiot*ratiot/(ratiol*ratiom);
		//dom.Particles[ip].bbeta = 0.1;
		//dom.Particles[ip].epsilon = 4.11e-20/(((ratiol/ratiot)*(ratiol/ratiot))*ratiom);
		//dom.Particles[ip].s = 1e-7/ratiol;
		//dom.Particles[ip].Lc = 6.28e-6/ratiol;
		//dom.Particles[ip].l = 3.04e-10/ratiol;
        //dom.Particles[ip].VdwCutoff = std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        dom.Particles[ip].D = 2;
		dom.Particles[ip].Ff =  0.0, -M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*my_dat.g,0.0, 0.0;
        //dom.Particles[ip].R = R;
        //dom.Particles[ip].M = M_PI*dom.Particles[ip].R*dom.Particles[ip].R*rhos;
        if(dom.Particles[ip].IsFree())
        {
            dom.Particles[ip].Rh = 0.8*R;
        }else{
            dom.Particles[ip].Rh = 0.8*RR;
        }  
    }
    double Tf = 100000;
    dom.IsF = true;    
    double dtout = 100;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolveIBM( Tf, dtout, "testclog", Setup, NULL);
    
}MECHSYS_CATCH
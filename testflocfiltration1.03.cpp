#include "./lbm/Domain.h"
////////////////////////////////////////////////////////////////////
//I have tested the max pdx = pdy = 60, and the pressure drop = 10 kPa. This combined parameters make the velcosity to be max to 0.46. So no much more than these parameters.
struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    double rhos;
    double rho1;
    double rho2;
    double sy;
    int    pnx;
    int    pny;
    double ratiol;
    double ratiot;
	double ratiom;
	size_t ny;
	double l;
	double lll;
	int    nfloc;
    int    fnx;
    int    shape;
    double RRR;
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
        if(dom.Particles[ip].X(1)>dat.ny-dat.lll-dat.R-2)
        {
            #pragma omp atomic
            ipp++;
        }
	}
    Vec3_t pos(0,0,0);
    Vec3_t v(0,0,0);
    Vec3_t w(0,0,0);
    int pnum = dom.Particles.size();
    if(dat.shape==3)
    {
        if(ipp<0)
        {
            for(int i=0; i<dat.fnx; ++i)
            {
                Vec3_t dxr(random(-0.3*dat.RRR,0.3*dat.RRR),random(0*dat.RRR,0*dat.RRR),0.0);
                int nnfloc = dat.nfloc;
                for(int k=-dat.nfloc; k<dat.nfloc+2; k=k+2)
                {
                    for (int kk=0; kk<nnfloc+1; ++kk)
                    {
                        double fx = 0.5*dat.lll + i*dat.lll + k*dat.R + kk*dat.R;
                        double fy =  dat.ny- dat.RRR - 2 - (std::sqrt(3)/3)*dat.nfloc*dat.R+ std::sqrt(3)*dat.R*kk;
                        pos = fx, fy , 0;
	                    dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, dat.rhos, dat.R, dom.dtdem));
                        pnum++;
                        dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
                        dom.Particles.back().Kn = 50.0;
                        dom.Particles.back().Gn = -0.3;
                        dom.Particles.back().Kt = 0.0;
                        dom.Particles.back().Mu = 0.4;
                        dom.Particles.back().Eta = 0.0;
                        dom.Particles.back().Beta = 0.0;
                        dom.Particles.back().A = 3e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
                        dom.Particles.back().kappa = 1e9*dat.ratiol;
                        dom.Particles.back().Z = 3.97e-11*dat.ratiot*dat.ratiot/(dat.ratiol*dat.ratiom);
    	                dom.Particles.back().bbeta = 0.1;
    	                dom.Particles.back().epsilon = 4.11e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
    	                dom.Particles.back().s = 1e-7/dat.ratiol;
    	                dom.Particles.back().Lc = 1.71e-7/dat.ratiol;
    	                dom.Particles.back().l = 3.04e-10/dat.ratiol;
                        dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
                        dom.Particles.back().D = 2;
                        dom.Particles.back().Rh = 0.80*dat.R;
                    }
                    nnfloc--; 
                } 
            }
        }
    }
    if(dat.shape==4)
    {
        if(ipp<0)
        {    
            for(int i=0; i<dat.fnx; ++i)
            {
                Vec3_t dxr(random(-0.45*dat.RRR,0.45*dat.RRR),random(0*dat.RRR,0*dat.RRR),0.0);
                for(int k=-dat.nfloc; k<dat.nfloc+2; k=k+2)
                {
		    		for(int kk=0; kk<dat.nfloc+1; ++kk)
                    {
                        double fx = 0.5*dat.lll+i*dat.lll + k*dat.R;
		    		    double fy = dat.ny - dat.RRR - 2 + (-dat.nfloc)*dat.R+2*dat.R*kk;
		    		    pos = fx, fy , 0;
		    		    dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, dat.rhos, dat.R, dom.dtdem));
		    		    pnum++;
                        dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
                        dom.Particles.back().Kn = 50.0;
                        dom.Particles.back().Gn = -0.3;
                        dom.Particles.back().Kt = 0.0;
                        dom.Particles.back().Mu = 0.4;
                        dom.Particles.back().Eta = 0.0;
                        dom.Particles.back().Beta = 0.0;
                        dom.Particles.back().A = 3e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
                        dom.Particles.back().kappa = 1e9*dat.ratiol;
                        dom.Particles.back().Z = 3.97e-11*dat.ratiot*dat.ratiot/(dat.ratiol*dat.ratiom);
    	                dom.Particles.back().bbeta = 0.1;
    	                dom.Particles.back().epsilon = 4.11e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
    	                dom.Particles.back().s = 1e-7/dat.ratiol;
    	                dom.Particles.back().Lc = 2.14e-6/dat.ratiol;
    	                dom.Particles.back().l = 3.04e-10/dat.ratiol;
                        dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
                        dom.Particles.back().D = 2;
                        dom.Particles.back().Rh = 0.80*dat.R;
                    }  
                }
            }  
        }
    }  
    if(dat.shape==6)
    {
        if(ipp<0)
        {
            double ffy = 0;
            double fy = dat.ny - dat.RRR - 2;
            for(int i=0; i<dat.fnx; ++i)
            {
                Vec3_t dxr(random(-0.45*dat.RRR,0.45*dat.RRR),random(0*dat.RRR,0*dat.RRR),0.0);
                for(int k=-dat.nfloc; k<(dat.nfloc+1); ++k)
                {
		    		double fx = 0.5*dat.lll+i*dat.lll+ k*2*dat.R;
                    pos = fx, fy , 0;
                    dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, dat.rhos, dat.R, dom.dtdem));
		    		pnum++;
                    dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
                    dom.Particles.back().Kn = 50.0;
                    dom.Particles.back().Gn = -0.3;
                    dom.Particles.back().Kt = 0.0;
                    dom.Particles.back().Mu = 0.4;
                    dom.Particles.back().Eta = 0.0;
                    dom.Particles.back().Beta = 0.0;
                    dom.Particles.back().A = 3e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
                    dom.Particles.back().kappa = 1e9*dat.ratiol;
                    dom.Particles.back().Z = 3.97e-11*dat.ratiot*dat.ratiot/(dat.ratiol*dat.ratiom);
    	            dom.Particles.back().bbeta = 0.1;
    	            dom.Particles.back().epsilon = 4.11e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
    	            dom.Particles.back().s = 1e-7/dat.ratiol;
    	            dom.Particles.back().Lc = 2.14e-6/dat.ratiol;
    	            dom.Particles.back().l = 3.04e-10/dat.ratiol;
                    dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
                    dom.Particles.back().D = 2;
                    dom.Particles.back().Rh = 0.80*dat.R;
		    		if(k == dat.nfloc)
		    		{
		    			continue;
		    		}
		    		if (k <= 0)
		    		{	
		    			for(int z =1; z<dat.nfloc+1; ++z)
		    			{
		    				fx = 0.5*dat.lll+i*dat.lll + k*2*dat.R + z*dat.R + dxr[0];
		    				ffy = fy + 2*dat.R*z*(std::sqrt(3))/2 ;
		    				pos = fx, ffy , 0;
		    				dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, dat.rhos, dat.R, dom.dtdem));
		    				pnum++;
                            dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
                            dom.Particles.back().Kn = 50.0;
                            dom.Particles.back().Gn = -0.3;
                            dom.Particles.back().Kt = 0.0;
                            dom.Particles.back().Mu = 0.4;
                            dom.Particles.back().Eta = 0.0;
                            dom.Particles.back().Beta = 0.0;
                            dom.Particles.back().A = 3e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
                            dom.Particles.back().kappa = 1e9*dat.ratiol;
                            dom.Particles.back().Z = 3.97e-11*dat.ratiot*dat.ratiot/(dat.ratiol*dat.ratiom);
    	                    dom.Particles.back().bbeta = 0.1;
    	                    dom.Particles.back().epsilon = 4.11e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
    	                    dom.Particles.back().s = 1e-7/dat.ratiol;
    	                    dom.Particles.back().Lc = 2.14e-6/dat.ratiol;
    	                    dom.Particles.back().l = 3.04e-10/dat.ratiol;
                            dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
                            dom.Particles.back().D = 2;
                            dom.Particles.back().Rh = 0.80*dat.R;
		    				ffy = fy - 2*dat.R*z*(std::sqrt(3))/2;
		    				pos = fx, ffy , 0;
		    				dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, dat.rhos, dat.R, dom.dtdem));
		    				pnum++;
                            dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
                            dom.Particles.back().Kn = 50.0;
                            dom.Particles.back().Gn = -0.3;
                            dom.Particles.back().Kt = 0.0;
                            dom.Particles.back().Mu = 0.4;
                            dom.Particles.back().Eta = 0.0;
                            dom.Particles.back().Beta = 0.0;
                            dom.Particles.back().A = 3e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
                            dom.Particles.back().kappa = 1e9*dat.ratiol;
                            dom.Particles.back().Z = 3.97e-11*dat.ratiot*dat.ratiot/(dat.ratiol*dat.ratiom);
    	                    dom.Particles.back().bbeta = 0.1;
    	                    dom.Particles.back().epsilon = 4.11e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
    	                    dom.Particles.back().s = 1e-7/dat.ratiol;
    	                    dom.Particles.back().Lc = 2.14e-6/dat.ratiol;
    	                    dom.Particles.back().l = 3.04e-10/dat.ratiol;
                            dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
                            dom.Particles.back().D = 2;
                            dom.Particles.back().Rh = 0.80*dat.R;
		    			}
		    		}
		    		else
		    		{
		    			for(int z =1; z<dat.nfloc+1-k; ++z)
		    			{
		    				fx = 0.5*dat.lll+i*dat.lll + k*2*dat.R + z*dat.R + dxr[0];
		    				ffy = fy + 2*dat.R*z*(std::sqrt(3))/2;
		    				pos = fx, ffy , 0;
		    				dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, dat.rhos, dat.R, dom.dtdem));
		    				pnum++;
                            dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
                            dom.Particles.back().Kn = 50.0;
                            dom.Particles.back().Gn = -0.3;
                            dom.Particles.back().Kt = 0.0;
                            dom.Particles.back().Mu = 0.4;
                            dom.Particles.back().Eta = 0.0;
                            dom.Particles.back().Beta = 0.0;
                            dom.Particles.back().A = 3e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
                            dom.Particles.back().kappa = 1e9*dat.ratiol;
                            dom.Particles.back().Z = 3.97e-11*dat.ratiot*dat.ratiot/(dat.ratiol*dat.ratiom);
    	                    dom.Particles.back().bbeta = 0.1;
    	                    dom.Particles.back().epsilon = 4.11e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
    	                    dom.Particles.back().s = 1e-7/dat.ratiol;
    	                    dom.Particles.back().Lc = 2.14e-6/dat.ratiol;
    	                    dom.Particles.back().l = 3.04e-10/dat.ratiol;
                            dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
                            dom.Particles.back().D = 2;
                            dom.Particles.back().Rh = 0.80*dat.R;
		    				ffy = fy - 2*dat.R*z*(std::sqrt(3))/2;
		    				pos = fx, ffy , 0;
		    				dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, dat.rhos, dat.R, dom.dtdem));
		    				pnum++;
                            dom.Particles.back().Ff =  0.0, -M_PI*dat.R*dat.R*(dat.rhos/1.0-1)*dat.g,0.0, 0.0;
                            dom.Particles.back().Kn = 50.0;
                            dom.Particles.back().Gn = -0.3;
                            dom.Particles.back().Kt = 0.0;
                            dom.Particles.back().Mu = 0.4;
                            dom.Particles.back().Eta = 0.0;
                            dom.Particles.back().Beta = 0.0;
                            dom.Particles.back().A = 3e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
                            dom.Particles.back().kappa = 1e9*dat.ratiol;
                            dom.Particles.back().Z = 3.97e-11*dat.ratiot*dat.ratiot/(dat.ratiol*dat.ratiom);
    	                    dom.Particles.back().bbeta = 0.1;
    	                    dom.Particles.back().epsilon = 4.11e-20/(((dat.ratiol/dat.ratiot)*(dat.ratiol/dat.ratiot))*dat.ratiom);
    	                    dom.Particles.back().s = 1e-7/dat.ratiol;
    	                    dom.Particles.back().Lc = 2.14e-6/dat.ratiol;
    	                    dom.Particles.back().l = 3.04e-10/dat.ratiol;
                            dom.Particles.back().VdwCutoff = std::sqrt(dom.Particles.back().A/(12.0*dom.Particles.back().Z*dom.Particles.back().kappa));
                            dom.Particles.back().D = 2;
                            dom.Particles.back().Rh = 0.80*dat.R;
		    			}
		    		}
                }
            } 
        }
    }
}

///////////////give the fluid an initialize speed,if give 0, not change these code until next comment/////////////////////////////////////////////////////
void Initial(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        // Vec3_t vtemp((double) dat.vb*iy/(ny-1), 0, 0);
        // Vec3_t vtemp((double) dat.vb, 0, 0);
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
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv) try
{
    std::srand((unsigned)time(NULL)); 
    size_t Nproc = 1;
    if    (argc>=2) Nproc = atoi(argv[1]); 
    printf("\n%s--- test n:  ---------------------------------------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // the particle   
    double nu = 0.10;
    int    Rn = 10;
    double R = Rn*1.0;                          //particle radius in the slurry
    int    pnx = 50;                            //small particle numbers in one row
    int    pny = 30;                            //the column number of small particles 
    double l  = 11;                             //gap between particle surfaces
    double rhos = 2.7;
    double Ga = 10.0;                           // Immersed boundary lattice boltzmann simulation of turbulent channel flows in the presence of spherical particles. Numerical Simulation of Two Spheres with Different Density Settling In Fluids through Lattice Boltzmann Method, 文章中定义了伽利略数的计算公式
 // the filter   
    double ratio = 2;
    double RR = ratio*R;                        //filter particle radius in the slurry
    int    Pnx = 1;                             //filter particle numbers in one row
    int    Pny = 1;                             //the column numbers of filter particles
    double pdx = 38.4;                            //gap between filter particle surface in x axis
    double pdy = 38.4;                           
 // the floc   
    int    fnx = 1;                             //floc number
    int    fny = 2;
    double lll = 126.5;                           //the gap between particle in the middle of the floc surfaces
    double nfloc = 2;                           //the number of the surface outside the particle in the middle of the floc
    double RRR = 0.0;                           //initialize the radius of the floc
    int    shape = 6;                           //shape =3 triangle; shape = 4 square; shape =6 hexagon
    if(shape==3)
    {
        RRR = (((nfloc*2)/3)*std::sqrt(3)+1)*R;
    }
    if(shape==4)
    {
        RRR = (nfloc*std::sqrt(2)+1)*R;
    }
    if(shape==6)
    {
        RRR = (nfloc*2+1)*R;
    }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      
 // the flow field
    size_t nx = std::max(std::ceil(Pnx*(RR*2+pdx)), std::ceil(fnx*lll));
	size_t ny = std::ceil((Pny-1)*(std::sqrt(3)*RR+pdy) + 40.0+2*RR + fny*lll)+1;
    size_t nz = 1;
    double rho = 1.0; 
    double dx = 1.0;
    double dt = 1.0;
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // reality convert to lattice
    double ratiol = 1e-6;                       //Lr/Ll
    double ratiot = 1e-7;                       //Tr/Tl 
	double ratiom = 1e-15;                      //Mr/Ml    	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                            
 // output the basic parameters in the test 
    std::cout<<"R   = "<<R<<std::endl;
    std::cout<<"RR  = "<<RR<<std::endl;
    std::cout<<"RRR = "<<RRR<<std::endl;
    std::cout<<"AOS/DDD = "<<(pdx/(RRR*2))<<std::endl;
    std::cout<<"AOS/D   = "<<(pdx/(R*2))<<std::endl;
    std::cout<<"nx1 = "<<std::ceil(Pnx*(RR*2+pdx))<<std::endl;
    std::cout<<"nx2 = "<<std::ceil(pnx*(R*2+l))<<std::endl;
    std::cout<<"nx3 = "<<std::ceil(fnx*lll)<<std::endl;
    std::cout<<"gap between floc surface = "<<std::ceil(lll-2*RRR)<<std::endl;   ///not precise
    if(shape == 3) 
    {
       
       std::cout<<"water content of the flocs = "<<((lll*lll-3.14*R*R*((nfloc+2)*(nfloc+1)/2))/(3.14*R*R*((nfloc+2)*(nfloc+1)/2)*2.7))<<std::endl; 
    }
    if(shape==4)
    {
        std::cout<<"water content of the flocs = "<<((lll*lll-3.14*R*R*(nfloc+1)*(nfloc+1))/((3.14*R*R*(nfloc+1)*(nfloc+1))*2.7))<<std::endl; 
    }
    if(shape==6)
    {
        std::cout<<"water content of the flocs = "<<((lll*lll-3.14*R*R*(3*(1+nfloc)*nfloc+1))/(3.14*R*R*(3*(1+nfloc)*nfloc+1)*2.7))<<std::endl;
    }
    std::cout<<"nx = "<<nx<<std::endl;
    std::cout<<"ny = "<<ny<<std::endl;
    std::cout<<"gy = "<<gy<<std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);         //去Domain.h文件看一段语句。
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = gy;
    my_dat.R = R;
    my_dat.rho1 = 1.03;
    my_dat.rho2 = 0.97;
    my_dat.pnx = pnx;
    my_dat.pny = pny;
	my_dat.ny = ny;
	my_dat.l = l;
    my_dat.rhos = rhos;
	my_dat.ratiol = ratiol;
	my_dat.ratiom = ratiom;
	my_dat.ratiot = ratiot;
	my_dat.lll = lll;
	my_dat.nfloc = nfloc;
    my_dat.rhos = rhos;
    my_dat.fnx = fnx;
    my_dat.shape = shape;
    my_dat.RRR = RRR;
    std::cout<<"pressure drop = "<<100*(my_dat.rho1-my_dat.rho2)/3<<"kPa"<<std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    dom.Nproc = Nproc;  
    dom.dtdem = 0.01*dt;
	Initial(dom,dom.UserData);
    bool iscontinue = false;
    Vec3_t v0(0.0,0.0,0.0);
    Vec3_t g0(0.0,0.0,0.0);
	//dom.InitialFromH5("testfloc_0043.h5",g0);
    //dom.Initial(rho,v0,g0);
    //bool iscontinue = true;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(!iscontinue)
    {	
		Vec3_t pos(0.0,0.0,0.0);
		Vec3_t pos1(0.0,0.0,0.0);
		Vec3_t dxp(0.0,0.0,0.0);
		Vec3_t v(0.0,0.0,0.0);
		Vec3_t w(0.0,0.0,0.0);
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    //filter particles and draw
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//floc position and draw
        double fy = 0;
        double fx = 0;
        // triangle floc
        if(shape == 3)
        {
            double sfy = (Pny-1)*(std::sqrt(3)*RR+pdy) + 2*RR + 0.5*lll + 40 + 1;
            //std::cout<<"sfy  = "<<sfy<<std::endl;
            for(int j=0; j<fny; ++j)
            {
                for(int i=0; i<fnx; ++i)
                {
                    Vec3_t dxr(random(-0.3*RRR,0.3*RRR),random(0*RRR,0*RRR),0.0);
                    int nnfloc = nfloc;
                    for(int k=-nfloc; k<nfloc+2; k=k+2)
                    {
                        for (int kk=0; kk<nnfloc+1; ++kk)
                        {
                            fx = 0.5*lll + i*lll + k*R + kk*R;
                            fy = sfy + j*lll - (std::sqrt(3)/3)*nfloc*R + std::sqrt(3)*R*kk;
                            //std::cout<<"fy  = "<<fy<<std::endl;
                            pos = fx, fy , 0;
	                        dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
                            pnum++;
                        }
                        nnfloc--; 
                    } 
                }
            }
        }
        //// square floc
        if(shape == 4)
        {
            double sfy = (Pny-1)*(std::sqrt(3)*RR+pdy) + 2*RR + 0.5*lll + 40 + 1;
            for(int j=0; j<fny; ++j)
            {
                for(int i=0; i<fnx; ++i)
                {
                    Vec3_t dxr(random(-0.45*RRR,0.45*RRR),random(0*RRR,0*RRR),0.0);
                    for(int k=-nfloc; k<nfloc+2; k=k+2)
                    {
		    			for(int kk=0; kk<nfloc+1; ++kk)
                        {
                            fx = 0.5*lll+i*lll + k*R;
		    			    fy = sfy + j*lll + (-nfloc)*R+2*R*kk;
		    			    pos = fx, fy , 0;
		    			    dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
		    			    pnum++;
                        }
                    }
                }
            }
        }
        //// hexagon floc
        if(shape == 6)
        {
            double ffy = 0;
            double sfy = (Pny-1)*(std::sqrt(3)*RR+pdy) + 2*RR + 0.5*lll + 40 + 1;
            for(int j=0; j<fny; ++j)
            {
                fy = sfy + j*lll;
                for(int i=0; i<fnx; ++i)
                {
                    Vec3_t dxr(random(-0.45*RRR,0.45*RRR),random(0*RRR,0*RRR),0.0);
                    for(int k=-nfloc; k<(nfloc+1); ++k)
                    {
		    			fx = 0.5*lll+i*lll+ k*2*R;
                        pos = fx, fy , 0;
                        dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
		    			pnum++;
		    			if(k == nfloc)
		    			{
		    				continue;
		    			}
		    			if (k <= 0)
		    			{	
		    				for(int z =1; z<nfloc+1; ++z)
		    				{
		    					fx = 0.5*lll+i*lll + k*2*R + z*R + dxr[0];
		    					ffy = sfy + j*lll + 2*R*z*(std::sqrt(3))/2 + dxr[1];
		    					pos = fx, ffy , 0;
		    					dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
		    					pnum++;
		    					ffy = sfy + j*lll - 2*R*z*(std::sqrt(3))/2 + dxr[1];
		    					pos = fx, ffy , 0;
		    					dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
		    					pnum++;
		    				}
		    			}
		    			else
		    			{
		    				for(int z =1; z<nfloc+1-k; ++z)
		    				{
		    					fx = 0.5*lll+i*lll + k*2*R + z*R + dxr[0];
		    					ffy = sfy + j*lll + 2*R*z*(std::sqrt(3))/2 + dxr[1];
		    					pos = fx, ffy , 0;
		    					dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
		    					pnum++;
		    					ffy = sfy + j*lll - 2*R*z*(std::sqrt(3))/2 + dxr[1];
		    					pos = fx, ffy , 0;
		    					dom.Particles.push_back(DEM::Disk(-pnum, pos, v, w, rhos, R, dom.dtdem));
		    					pnum++;
		    				}
		    			}
                    }
                }
            }
        }
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    //the basic parameters for filter particles and slurry particles
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Kn = 50.0;
        dom.Particles[ip].Gn = -0.3;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.4;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].A = 3e-20/(((ratiol/ratiot)*(ratiol/ratiot))*ratiom);
        dom.Particles[ip].kappa = 1e9*ratiol;
        dom.Particles[ip].Z = 3.97e-11*ratiot*ratiot/(ratiol*ratiom);
		dom.Particles[ip].bbeta = 0.1;
		dom.Particles[ip].epsilon = 4.11e-20/(((ratiol/ratiot)*(ratiol/ratiot))*ratiom);
		dom.Particles[ip].s = 1e-7/ratiol;
		dom.Particles[ip].Lc = 1.71e-7/ratiol;
		dom.Particles[ip].l = 3.04e-10/ratiol;
        dom.Particles[ip].VdwCutoff = std::sqrt(dom.Particles[ip].A/(12.0*dom.Particles[ip].Z*dom.Particles[ip].kappa));
        dom.Particles[ip].D = 2;
		dom.Particles[ip].Ff =  0.0, -M_PI*dom.Particles[ip].R*dom.Particles[ip].R*(rhos/rho-1)*my_dat.g,0.0, 0.0;
        if(dom.Particles[ip].IsFree())
        {
            dom.Particles[ip].Rh = 0.80*R;
        }
        else
        {
            dom.Particles[ip].Rh = 0.80*RR;
        }
    }
	std::cout<<"VdwCutoff "<<dom.Particles.back().VdwCutoff<<std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double Tf = 100000;
	std::cout<<"9 Tf="<<Tf<<std::endl;
    dom.IsF = true;    
    double dtout = 100;
    // periodic in x
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    // periodic in y
    // dom.Box = 0.0, ny-1, 0.0;
    // dom.modexy = 1;
    //solving
    dom.SolveIBM( Tf, dtout, "D4_", Setup, NULL);
}MECHSYS_CATCH
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <objective_function.hpp>
#include <EMC.h>
#include <mpi.h>

int nfields = 3;
int ndim = 2;
int nx = 50;
int ny = 50;
double xmax = 10.0;
double ymax = 10.0;
double dx,dy;
int testsum = 0;

double B = 1.0;

double kappa = 1.0;
double m_pi = sqrt(0.1/2.0);
double charge_chem = 2.0;

int imap(int n,int i,int j)
{
    return nx*ny*n + ny*i + j;
}

int imapsmall(int n,int i,int j)
{
    return (nx-4)*(ny-4)*n + (ny-4)*i + j;
}

double profile(double r){
    return M_PI*exp(-0.06*pow(r,2));
}

double chargedensity(int i, int j, const  Eigen::ArrayXd & f){
    double chargeden;
    double df[ndim][nfields];

    for(int n=0;n<nfields;n++) {
        df[0][n] = (-f(imap(n,i+2,j)) + 8.0 * f(imap(n,i+1,j)) - 8.0 * f(imap(n,i-1,j)) +
                    f(imap(n,i-2,j)))/(12.0*dx);
        df[1][n] = (-f(imap(n,i,j+2)) + 8.0 * f(imap(n,i,j+1)) - 8.0 * f(imap(n,i,j-1)) +
                    f(imap(n,i,j-2)))/(12.0*dy);
    }


    chargeden = - (1.0/(4.0*M_PI))*f(imap(0,i,j))*(df[0][1]*df[1][2] - df[0][2]*df[1][1]);
    chargeden = chargeden - (1.0/(4.0*M_PI))*f(imap(1,i,j))*(-1.0)*(df[0][0]*df[1][2] - df[0][2]*df[1][0]);
    chargeden = chargeden - (1.0/(4.0*M_PI))*f(imap(2,i,j))*(df[0][0]*df[1][1] - df[0][1]*df[1][0]);

    return chargeden;
}

double energydensity(int i,int j, const Eigen::ArrayXd  & f){
    double enden = 0.0;
    double df[ndim][nfields];
    double no_good;
    no_good = 0.0;

    for(int n=0;n<nfields;n++) {
        df[0][n] = (-f(imap(n,i+2,j)) + 8.0 * f(imap(n,i+1,j)) - 8.0 * f(imap(n,i-1,j)) +
                    f(imap(n,i-2,j)))/(12.0*dx);
        df[1][n] = (-f(imap(n,i,j+2)) + 8.0 * f(imap(n,i,j+1)) - 8.0 * f(imap(n,i,j-1)) +
                    f(imap(n,i,j-2)))/(12.0*dy);
    }

    for(int n=0;n<nfields;n++) {
        for(int mu=0;mu<ndim;mu++){
            enden = enden + 0.25*df[mu][n]*df[mu][n];
            for(int l=0;l<nfields;l++){
                for(int nu=0;nu<ndim;nu++){
                    enden = enden + (pow(kappa,2)/4.0)*(df[mu][n]*df[mu][n]*df[nu][l]*df[nu][l] - df[mu][n]*df[nu][n]*df[mu][l]*df[nu][l]);
                }

            }
        }
    }

    return enden + pow(m_pi,2)*( pow(f(imap(0,i,j)),2) + pow(f(imap(1,i,j)),2) );
}

double sum_energy (objective_function &obj_fun, Eigen::ArrayXd & in1){
    testsum = testsum+1;
    auto in = TensorMap<Tensor<double,3>>(in1.data(),nfields,nx-4,ny-4);  //(3 is rank, and 3,3,3 the sizes of the tensor)
    double energy = 0.0;
    double charge = 0.0;
    Eigen::ArrayXd  f(nfields*nx*ny);
    f.fill(1);
    for(int i=0;i<nx;i++) {
        for (int j=0;j<ny;j++) {
            if (i < 2 || i >= nx - 2 || j < 2 || j >= ny - 2) {
                //set vacuum
                f(imap(0, i, j)) = 0.0;
                f(imap(1, i, j)) = 0.0;
                f(imap(2, i, j)) = 1.0;
            } else {
                //renormalise
                double r = sqrt(pow(in(0, i-2, j-2), 2) + pow(in(1, i-2, j-2), 2) + pow(in(2, i-2, j-2), 2));
                if(r > 0) {
                    f(imap(0, i, j)) = in(0, i-2, j-2) / r;
                    f(imap(1, i, j)) = in(1, i-2, j-2) / r;
                    f(imap(2, i, j)) = in(2, i-2, j-2) / r;
                }else{cout << "paramterss are all zero!!!\n";}
            }
        }
    }

    //want to loop through the field theory and find some density based on derivatives of the various fields
    for(int i=2;i<nx-2;i++)
    {
        for(int j=2;j<ny-2;j++)
        {
            charge += dx*dy*chargedensity(i,j,f);
            energy += dx*dy*energydensity(i,j,f);
        }
    }
    //cout << "returning en\n";
    return (energy/(2.0*M_PI)) + charge_chem*fabs(B - charge);
   // return sqrt(inputParameters.cwiseAbs2().sum())+10.0*abs(sin(5.0*sqrt(abs(inputParameters.cwiseAbs2().sum())))) ;
}

double sum_energy_final (Eigen::Tensor<double,3> & in){
    //auto in = TensorMap<Tensor<double,3>>(in1.data(),nfields,nx-4,ny-4);  //(3 is rank, and 3,3,3 the sizes of the tensor)
    double energy = 0.0;
    double charge = 0.0;
    Eigen::ArrayXd  f(nfields*nx*ny);
    f.fill(1);
    for(int i=0;i<nx;i++) {
        for (int j=0;j<ny;j++) {
            if (i < 2 || i >= nx - 2 || j < 2 || j >= ny - 2) {
                //set vacuum
                f(imap(0, i, j)) = 0.0;
                f(imap(1, i, j)) = 0.0;
                f(imap(2, i, j)) = 1.0;
            } else {
                //renormalise
                double r = sqrt(pow(in(0, i-2, j-2), 2) + pow(in(1, i-2, j-2), 2) + pow(in(2, i-2, j-2), 2));
                if(r > 0) {
                    f(imap(0, i, j)) = in(0, i-2, j-2) / r;
                    f(imap(1, i, j)) = in(1, i-2, j-2) / r;
                    f(imap(2, i, j)) = in(2, i-2, j-2) / r;
                }else{cout << "paramters are all zero!!!\n";}
            }
        }
    }

    //want to loop through the field theory and find some density based on derivatives of the various fields
    for(int i=2;i<nx-2;i++)
    {
        for(int j=2;j<ny-2;j++)
        {
            charge += dx*dy*chargedensity(i,j,f);
            energy += dx*dy*energydensity(i,j,f);
        }
    }
    //cout << "returning en\n";
    return (energy/(2.0*M_PI));
    // return sqrt(inputParameters.cwiseAbs2().sum())+10.0*abs(sin(5.0*sqrt(abs(inputParameters.cwiseAbs2().sum())))) ;
}

double sum_charge_final (Eigen::Tensor<double,3> & in){
    //auto in = TensorMap<Tensor<double,3>>(in1.data(),nfields,nx-4,ny-4);  //(3 is rank, and 3,3,3 the sizes of the tensor)
    double energy = 0.0;
    double charge = 0.0;
    Eigen::ArrayXd  f(nfields*nx*ny);
    f.fill(0);
    for(int i=0;i<nx;i++) {
        for (int j=0;j<ny;j++) {
            if (i < 2 || i >= nx - 2 || j < 2 || j >= ny - 2) {
                //set vacuum
                f(imap(0, i, j)) = 0.0;
                f(imap(1, i, j)) = 0.0;
                f(imap(2, i, j)) = 1.0;
            } else {
                //renormalise
                double r = sqrt(pow(in(0, i-2, j-2), 2) + pow(in(1, i-2, j-2), 2) + pow(in(2, i-2, j-2), 2));
                if(r > 0) {
                    f(imap(0, i, j)) = in(0, i-2, j-2) / r;
                    f(imap(1, i, j)) = in(1, i-2, j-2) / r;
                    f(imap(2, i, j)) = in(2, i-2, j-2) / r;
                }else{cout << "paramters are all zero!!!\n";}
            }
        }
    }

    //want to loop through the field theory and find some density based on derivatives of the various fields
    for(int i=2;i<nx-2;i++)
    {
        for(int j=2;j<ny-2;j++)
        {
            charge += dx*dy*chargedensity(i,j,f);
            energy += dx*dy*energydensity(i,j,f);
        }
    }
    //cout << "returning en\n";
    return charge;
    // return sqrt(inputParameters.cwiseAbs2().sum())+10.0*abs(sin(5.0*sqrt(abs(inputParameters.cwiseAbs2().sum())))) ;
}

int main() {
    dx = 2.0*xmax/nx;
    dy = 2.0*ymax/ny;
    Eigen::Tensor<double,3> lower_bound(nfields,nx-4,ny-4);
    Eigen::Tensor<double,3> upper_bound(nfields,nx-4,ny-4);
    Eigen::Tensor<double,3> initialconditions(nfields,nx-4,ny-4);
    for(int i=0;i<nx-4;i++) {
        for (int j=0;j<ny-4;j++) {
            double r = sqrt(pow(2.0*(i+2.0)*dx-xmax,2)+pow(2.0*(j+2.0)*dy-ymax,2));
            double theta = atan2(2.0*(j+2.0)*dy-ymax,2.0*(i+2.0)*dx-xmax);
            initialconditions(0, i, j) = sin(profile(r)) * cos(theta);
            initialconditions(1, i, j) = sin(profile(r)) * sin(theta);
            initialconditions(2, i, j) = cos(profile(r));
            //cout << "i=" << i << " j=" << j << " x=" << 2.0*(i+2.0)*dx - xmax << " y= " << 2.0*(j+2.0)*dy-ymax << " r=" << r << " theta=" << theta << " f0=" << f(imap(0,i,j))<< " f1=" << f(imap(1,i,j))<< " f2=" << f(imap(2,i,j))<<"\n";
        }
    }
    lower_bound.setConstant(-1.0);
    upper_bound.setConstant(1.0);
    objective_function obj_fun(sum_energy,lower_bound,upper_bound,1e-8,initialconditions);
//    objective_function obj_fun(sum_energy,lower_bound,upper_bound,1e-6,ArrayXd::Zero(1));
    cout << "initial energy = " << sum_energy_final(initialconditions) << " with charge = " << sum_charge_final(initialconditions)<<"\n";
    //obj_fun.threads=6;
    minimize(obj_fun);
    //std::cout << obj_fun.optimum << std::endl;
    //cout << "final energy = " << sum_energy_final(obj_fun.optimum) << " with charge = " << sum_charge_final(obj_fun.optimum)<<"\n";
    return 0;
}
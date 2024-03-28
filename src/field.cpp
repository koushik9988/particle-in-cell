#include "field.h"


//void FieldSolve::SolvePotDirect(double *x, double *rho)
void FieldSolve::SolvePotDirect()
{
    /* Set coefficients, precompute them*/
    double *x = domain.phi;
    double *rho = domain.rho;
    int ni = domain.ni;
    double dx2 = domain.dx*domain.dx;
    double *a = new double[ni];
    double *b = new double[ni];
    double *c = new double[ni];

    memset(a, 0, ni * sizeof(double));
    memset(b, 0, ni * sizeof(double));
    memset(c, 0, ni * sizeof(double));

    /*Centtral difference on internal nodes*/
    for(int i=1; i<ni-1; i++)
    {
        a[i] = 1; b[i] = -2; c[i] = 1;
    }

    /*Apply dirichlet boundary conditions on boundaries*/
    a[0]= 0; b[0]=1; c[0]=0;
    a[ni-1]=0; b[ni-1]=1; c[ni-1]=0;

    /*multiply R.H.S.*/
    for (int i=1; i<ni-1; i++)
        x[i]=-rho[i]*dx2;

    //boundary conditions
    //x[0] = 0;
    //x[ni-1] = 0;
    x[0] = x[ni-1];

    /*Modify the coefficients*/
    c[0] /=b[0];
    x[0] /=b[0];

    for(int i=1; i<ni; i++)
    {
        double id = (b[i]-c[i-1]*a[i]);
        c[i] /= id;
        x[i] = (x[i]-x[i-1]*a[i])/id;
    }

    /* Now back substitute */
    for(int i=ni-2; i>=0; i--)
        x[i] = x[i] - c[i]*x[i+1];

    //delete allocated memory
    delete [] a;
    delete [] b;
    delete [] c;
    //return x;
}

bool FieldSolve::SolvePotIter()
{
    double L2;
    double *phi = domain.phi;
    double *rho = domain.rho;
    int ni = domain.ni;
    double dx2 = domain.dx*domain.dx;

    // Initialize boundaries
    phi[0]=phi[ni-1] = 0;

    // Main Solver
    for(int it=0; it<40000; it++)
    {
        for(int i=1; i< ni-1; i++)
        {
            double g = 0.5*(phi[i-1] + phi[i+1] + dx2*rho[i]);
            phi[i]=phi[i] + 1.4*(g-phi[i]);
        }
        // Check for convergence
        if(it%25==0)
        {
            double sum = 0;
            for(int i=1; i< ni-1; i++)
            {
                double R = - rho[i] - (phi[i-1]-2*phi[i]+phi[i+1])/dx2;
                sum += R*R;
            }
            L2 = sqrt(sum)/domain.ni;
            if(L2<1e-4){return true;}

        }
        //printf("GS-Converged! L2=%g\n",L2);
    }
    printf("Gauss-Siedel solver failed to converge, L2=%g\n",L2);
    return false;
}

//void FieldSolve::CalculateEfield(double *phi, double *ef)
void FieldSolve::CalculateEfield()
{
    double *phi = domain.phi;
    double *ef = domain.ef;
    //std::string bc  = domain.bc;
    /*Apply central difference to the inner nodes*/
    for(int i=1; i<domain.ni-1; i++)
    {
         ef[i] = -(phi[i+1]-phi[i-1])/(2*domain.dx);
    }
       
    /*for continous bounndary
    the point 0 and ni-1 is same */
    if(domain.bc =="pbc")
    {
        ef[0] = -(phi[1]-phi[domain.ni-2])/(2*domain.dx);
        ef[domain.ni-1] = ef[0];
    }
    else if(domain.bc == "open")
    {
        ef[0] = -(phi[1]-phi[0])/(domain.dx);
        ef[domain.ni-1] = -(phi[domain.ni-1]-phi[domain.ni-2])/(domain.dx);
    }
    
    
}
#include "field.h"

//mo problem here(numerical error problem )
//void FieldSolve::SolvePotDirect(double *x, double *rho)
void FieldSolve::Direct(int ts)
{
    /* Set coefficients, precompute them*/
    int ni = domain.ni;
    double dx2 = domain.dx*domain.dx;

    vec<double> rhov(ni);

    //boundary condition

    
    double frequency = 13.56e6;
    double omega = 2.0 * Const::PI * frequency;


    double volatge = 100.0;  // Amplitude of the electric field pulse
    
    //double t_start = 10;  // Start time of the pulse 
    //double t_end = 5000;    // End time of the pulse
    //double current_time = ts * domain.DT;  // Current time in simulation

    if(domain.bc == "pbc")
    {  
        /*

        if (current_time >= t_start && current_time <= t_end)
        {
            // Apply pulse at the left boundary for the duration of t_start to t_end
            rhov(0) = 0*pulse_amplitude * sin(omega*ts*domain.DT) * (Const::eV / (Const::K_b*Const::EV_to_K));
            rhov(ni-1) = 0*pulse_amplitude * cos(omega*ts*domain.DT) * (Const::eV / (Const::K_b*Const::EV_to_K));
        }
        else
        {
            rhov(0) = 0;
            rhov(ni-1) = 0;
        }
        
        */
        //domain.vL =  volatge*sin(omega*ts*domain.DT)*(Const::eV/(Const::K_b*Const::EV_to_K));
        rhov(0) = 0;//domain.vL;// 10*sin(omega*ts*domain.DT*1)*(Const::eV/(Const::K_b*Const::EV_to_K));
        rhov(ni-1) = 0;//100*cos(omega*ts*domain.DT*10);   
    }
    if(domain.bc == "open")
    {
        //rhov(0) =  10*sin(omega*ts*domain.DT*10)*(Const::eV/(Const::K_b*Const::EV_to_K));//domain.vL;
        //domain.vL =  volatge*sin(omega*ts*domain.DT)*(Const::eV/(Const::K_b*Const::EV_to_K));
        //domain.vR = 10000*(Const::eV/(Const::K_b*Const::EV_to_K));
        rhov(0) = domain.vL; 
        rhov(ni-1) = domain.vR;//domain.vR;  
    }

                 
    for( int i = 1 ; i < ni-1 ; i++)
    {
        rhov(i) = -domain.rho(i)*dx2*((domain.L*domain.L)/(domain.LDe*domain.LDe));
        //rhov(i) += 0.1*sin((2*3.14*i)/domain.xL)*((domain.L*domain.L)/(domain.LD*domain.LD));
    }

    vec<double> sol = direct(rhov,ni);

    domain.phi = 0;
    domain.phi = sol;
}


void FieldSolve::GaussElim()
{
    double L2;
    int n = domain.ni;
    double dx2 = domain.dx*domain.dx;

    vec<double> b(n);
    Matrix<double> A(n,n);

   
    A(0, 0) = 1.0; //
    b(0) = 0;   //left fixed potential  bc
    
    A(n - 1, n - 1) = 1.0; 
    b(n - 1) = 0; // Right boundary (fixed potential)

    // Populate internal nodes
    for (int i = 1; i < n-1; i++) 
    {
        A(i, i - 1) = 1.0; // A(n, n-1)
        A(i, i) = -2.0;   // A(n, n)
        A(i, i + 1) = 1.0; // A(n, n+1)
        b(i) = -domain.rho(i)*dx2*((domain.L*domain.L)/(domain.LDe*domain.LDe));
    }
    
    vec<double> sol = gausselimination(A,b);

    domain.phi = 0;
    
    domain.phi = sol;
}


void FieldSolve::pcgsolver()
{
    double L2;
    int n = domain.ni;
    double dx2 = domain.dx*domain.dx;

    vec<double> b(n);
    vec<double> x(n);
    Matrix<double> A(n,n);

    std::vector<bool> fixed_node(n, false); 
    // Set boundary conditions
   
    A(0, 0) = 1.0; //
    b(0) = 0;   //left fixed potential  bc
    
    A(n - 1, n - 1) = 1.0; 
    b(n - 1) = 0; // Right boundary (fixed potential)

    // Populate internal nodes
    for (int i = 1; i < n-1; i++) 
    {
        A(i, i - 1) = 1.0; // A(n, n-1)
        A(i, i) = -2.0;   // A(n, n)
        A(i, i + 1) = 1.0; // A(n, n+1)
        b(i) = -domain.rho(i)*dx2*((domain.L*domain.L)/(domain.LDe*domain.LDe));
        x(i) = domain.phi(i);
    }

    Matrix<double> A0 = slice(A,1,n-1,1,n-1);
    vec<double> b0 = slice(b,1,n-1);
    vec<double> x0 = slice(x,1,n-1);

    //symmetrycheck(A0);

    vec<double> x_sol = cg(A0,x0,b0,domain.max_iteration, domain.tolerance);
   
    vec<double> v_norm = A0*x_sol - b0;

    //reset phi
    domain.phi = 0;

    domain.phi(0) = b(0);
    domain.phi(n-1) = b(n-1);
    
    domain.phi = x_sol;

    for(int i = 1; i< n-1; i++)
    {
        domain.phi(i) = x_sol(i-1);
    }
}


void FieldSolve::cgsolver()
{
    double L2;

    int n = domain.ni;
    double dx2 = domain.dx*domain.dx;

    vec<double> b(n);
    vec<double> x(n);
    Matrix<double> A(n,n);

    std::vector<bool> fixed_node(n, false); 
    // Set boundary conditions
   
    A(0, 0) = 1.0; //
    b(0) = 0;   //left fixed potential  bc
    
    A(n - 1, n - 1) = 1.0; 
    b(n - 1) = 0; // Right boundary (fixed potential)

    // Populate internal nodes
    for (int i = 1; i < n-1; i++) 
    {
        A(i, i - 1) = 1.0; // A(n, n-1)
        A(i, i) = -2.0;   // A(n, n)
        A(i, i + 1) = 1.0; // A(n, n+1)
        b(i) = -domain.rho(i)*dx2*((domain.L*domain.L)/(domain.LDe*domain.LDe));;
        x(i) = domain.phi(i);
    }

    Matrix<double> A0 = slice(A,1,n-1,1,n-1);
    vec<double> b0 = slice(b,1,n-1);
    vec<double> x0 = slice(x,1,n-1);

    //symmetrycheck(A0);

    vec<double> x_sol = cg(A0,x0,b0,domain.max_iteration, domain.tolerance);
   
    vec<double> v_norm = A0*x_sol - b0;

    //display::print(x_sol.getSize());

    domain.norm = v_norm.norm();

    domain.phi(0) = b(0);
    domain.phi(n-1) = b(n-1);
    
    for(int i = 1; i< n-1; i++)
    {
        domain.phi(i) = x_sol(i-1);
    }
}

//void FieldSolve::CalculateEfield(double *phi, double *ef)
void FieldSolve::CalculateEfield()
{
    //std::string bc  = domain.bc;
    /*Apply central difference to the inner nodes*/
    domain.ef = 0;
    for(int i=1; i<domain.ni-1; i++)
    {
        domain.ef(i) = -((domain.LDe*domain.LDe)/(domain.L*domain.L))*(domain.phi(i+1)-domain.phi(i-1))/(2*domain.dx);
    }
       
    /*for continous bounndary
    the point 0 and ni-1 is same */
    if(domain.bc =="pbc")
    {
        domain.ef(0) = -((domain.LDe*domain.LDe)/(domain.L*domain.L))*(domain.phi(1)-domain.phi(domain.ni-2))/(2*domain.dx);
        domain.ef(domain.ni-1) = domain.ef(0);
    }
    else if(domain.bc == "open" || domain.bc == "rbc")
    {
        domain.ef(0) = -((domain.LDe*domain.LDe)/(domain.L*domain.L))*(domain.phi(1)- domain.phi(0))/(domain.dx);
        domain.ef(domain.ni-1) = -((domain.LDe*domain.LDe)/(domain.L*domain.L))*(domain.phi(domain.ni-1)-domain.phi(domain.ni-2))/(domain.dx);
    }
}


void FieldSolve::PotentialSolver(int ts)
{
    std::string solverTypeStr = domain.SolverType;
    SolverType solverType;

    if (solverTypeStr == "direct")
        solverType = SolverType::DIRECT;
    else if (solverTypeStr == "pcg")
        solverType = SolverType::PCG;
    else if (solverTypeStr == "cg")
        solverType = SolverType::CG;
    else if (solverTypeStr == "gs")
        solverType = SolverType::GS;
    else
        solverType = SolverType::DIRECT; // Default to DIRECT

    switch (solverType)
    {
    case SolverType::DIRECT:
        Direct(ts);
        break;
    case SolverType::PCG:
        pcgsolver();
        break;
    case SolverType::CG:
        cgsolver();
        break;
    case SolverType::GS:
        GaussElim();
        break;
    default:
        Direct(ts);
        break;
    }
}
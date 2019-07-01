
    #include <iostream>
    #include </usr/local/include/gsl/gsl_sf_exp.h>
    #include </usr/local/include/gsl/gsl_errno.h>
    #include </usr/local/include/gsl/gsl_sf_coulomb.h>
    #include <stdlib.h>
    #include "sfactor3Hedp.h"
    const double cs_pi = 0.31415926535897932384626433832795028e+01;

    using std::vector;
sfactor3Hedp::sfactor3Hedp() : ArrayFunction("sfactor3Hedp", 8)
    {
    }
    void sfactor3Hedp::evaluate (double *value,
    std::vector<double const *> const &args,
    std::vector<std::vector<unsigned int> > const &dims) const
    {
    double E  = args[0][0];		//energy at which S-factor is to be calculated
    double e1 = args[1][0];		//eigenenergy
    double ex = args[2][0];     //location of sigma or S-factor peak
    double gi = args[3][0];		//reduced width of incoming channel
    double gf = args[4][0];		//reduced width of outgoing channel
    double ri = args[5][0];		//radius of incoming channel
    double rf = args[6][0];		//radius of of outgoing channel
    double ue = args[7][0];		//electron screening potential


    double m1_i = 3.01493, m2_i = 2.01355;	//masses (amu) of t and d
    double m1_f = 4.00151, m2_f = 1.007277;	//masses (amu) of n and 4He
    double z1_i = 2, z2_i = 1;				//charges of t and d
    double z1_f = 2, z2_f = 1;				//charges of n and 4He
    double jt=0.5, jp=1.0, jr=1.5;			//spins of target, projectile, resonance
    double Q = 18.353053;					//reaction Q-value (MeV)
    int    la = 0, lb = 2;					//orbital angular momenta of d and n


    double mue_i=(m1_i*m2_i)/(m1_i+m2_i);
    double mue_f=(m1_f*m2_f)/(m1_f+m2_f);
    double pek = 6.56618216e-1/mue_i;
    double omega = (2*jr+1)/( (2*jt+1)*(2*jp+1) );
    double s1, s2;
    double F,FP, G,GP;


    double etpe_i=exp( 0.98951013e0*z1_i*z2_i*sqrt(mue_i/E) );

    double p_i, px_i, s_i, b_i, Ga;
    double p_f, px_f, s_f, b_f, Gb;
    double tapp;
     /* INCOMING CHANNEL
     /* penetration and shift factors
    PenFactor(E, la, ri, mue_i, z1_i*z2_i, p_i, s_i);


     /* calculate boundary condition parameter
    PenFactor(ex, la, ri, mue_i, z1_i*z2_i, px_i, b_i);

    Ga = 2*gi*p_i; // Rafa

    /* OUTGOING CHANNEL
    PenFactor(E+Q, lb, rf, mue_f, z1_f*z2_f, p_f, s_f);


    /* calculate boundary condition parameter
    PenFactor(ex+Q, lb, rf, mue_f, z1_f*z2_f, px_f, b_f);


    Gb = 2*gf*p_f; // Rafa
     /* boundary condition: Sc(ex), where ex is near peak of observed sigma or S
    tapp=(s_i-b_i)*gi+(s_f-b_f)*gf;        // level shift
     /*     ---------------------------------------------------
     /*     PUTTING EVERYTHING TOGETHER
     /*     ----------------------------------------------------
    s1=pek*etpe_i*omega*Ga*Gb;
    s2=( pow(e1-E-tapp,2) )+0.25e0*( pow(Ga+Gb,2) );

     /* we output the S-factor by storing it in value[0]
    value[0] = exp( 0.5*0.98951013e0*z1_i*z2_i*sqrt(mue_i)*ue*pow(E,-1.5) )*s1/s2;
    }

     /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /* Function to calcuate the penetration factor
     /* INPUT: E, L, Mass0, Mass1, Charge0, Charge1
     /* OUTPUT: P, S
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void sfactor3Hedp::PenFactor(
    const double E,
    const double L,
    const double R,
    const double mue,
    const double qQ,
    double& P,
    double& S) const {

     /* Turn off the GSL error handler, which aborts the program
     /* if G or F go out of range.
    gsl_set_error_handler_off();

    gsl_sf_result F,Fp,G,Gp;

    double exp_F,exp_G;

    double eta = 0.1574854*qQ*sqrt(mue/E);
    double rho = 0.218735*R*sqrt(mue*E);

     /*this needs eta, rho, L
     /*everything above this is just to calculate eta, rho
    int status = gsl_sf_coulomb_wave_FG_e (eta, rho, L, 0, &F, &Fp, &G,
    &Gp, &exp_F, &exp_G);
    if(status){
      if(status == GSL_EOVRFLW){
        /*ErrorFlag = true;
        /*PenZeroCount++;
        P = 0.0;
        S = 0.0;
        return;
      }
      else {
        exit(1);
        std::cout << "
ERROR: Something went wrong in coulomb wavefunction!" << "
	GSL Error: " << gsl_strerror (status)<< std::endl;
        std::cout << "The Energy was " << E*1e3 << " keV." << std::endl;
        abort();
      }
    }
    gsl_set_error_handler (NULL);

     /* Just in case there is an overflow, multiply by exponential
     /* (See GSL documentation for more info)
    double F_l = F.val*exp(exp_F);
    double G_l = G.val*exp(exp_G);

    P = rho/( pow(F_l,2) + pow(G_l,2) );
    S = rho*( F_l * Fp.val + G_l * Gp.val)/(pow(F_l,2) + pow(G_l,2));
    return;
    }

/**
  * Checks whether dimensions of the function parameters are correct.
*
  * @param dims Vector of length npar denoting the dimensions of
* the parameters, with any redundant dimensions dropped.
*/
  bool sfactor3Hedp::checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) const {
     /*the first argument should be a vector
     /* the last three arguments should be scalars
    return isScalar(dims[0]) && isScalar(dims[1]) && isScalar(dims[2]) && isScalar(dims[3])
    && isScalar(dims[4]) && isScalar(dims[5]) && isScalar(dims[6])
    && isScalar(dims[7]);
  }
/**
  * Checks whether the parameter values lie in the domain of the
* function. The default implementation returns true.
*/
  bool sfactor3Hedp::checkParameterValue(std::vector<double const *> const &args,
                                         std::vector<std::vector<unsigned int> > const &dims) const{
                                           // TODO: should any parameters be eg strictly positive?
                                           return true;
                                         }

/**
  * Calculates what the dimension of the return value should be,
* based on the arguments.
*
  * @param dims Vector of Indices denoting the dimensions of the
* parameters. This vector must return true when passed to
* checkParameterDim.
*
  * @param values Vector of pointers to parameter values.
*/
  std::vector<unsigned int> sfactor3Hedp::dim(std::vector <std::vector<unsigned int> > const &dims,
                                              std::vector <double const *> const &values) const {
                                                // the size of the table that the fortran code calculates is length of E (input)
                                                vector<unsigned int> ans(1);

                                                ans[0] = 1;
                                                return ans;
                                              }

    
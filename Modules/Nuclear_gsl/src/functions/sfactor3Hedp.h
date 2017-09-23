#ifndef SFACTOR3HEDP_H_
#define SFACTOR3HEDP_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace nuclear {

/**
 * @short Astrophysical S Factors
 * factor3Hedp returns the astrophysical S-factor as a function of energy.  It
 * returns a 1000x2 table where column 1 is energy and column 2 is the S-factor.
 * It requires a fortran executable to perform the calculations.
 * <pre>
 * table = factor3Hedp(e1, gi, gf)
 * </pre>
 */

    class sfactor3Hedp : public ArrayFunction
    {
    public:
	sfactor3Hedp();
void evaluate(double *x, std::vector<double const *> const &args,
                    std::vector<std::vector<unsigned int> > const &dims) const;
                    
    void coul(int, double, double, double&, double&) const;
    
   /*void PenFactor(const double E, const double L, const double R,
   				const double mue,
				const double qQ,
				double& P, double& S) const;*/
    
                
	bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) const;
	
	std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims,
    	
    	std::vector<double const *> const &values) const;
    	
    	
	bool checkParameterValue(std::vector<double const *> const &args,
    	std::vector<std::vector<unsigned int> > const &dims) const;
}; 
    
    
}}//end namespaces

#endif /* SFACTOR3HEDP_H_ */

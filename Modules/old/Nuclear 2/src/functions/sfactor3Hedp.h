#ifndef SFACTOR3HEDP_H_
#define SFACTOR3HEDP_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace nuclear {

/**
 * @short Astrophysical S Factors
 * sfactor3Hedp_5p returns the astrophysical S-factor as a function of energy.
 * <pre>
 * obsy1 = sfactor3Hedp(obsx1, e1, gi, gf, ri, rf, ue)
 * </pre>
 */

class sfactor3Hedp : public ArrayFunction {
    public:
	sfactor3Hedp();
	
    void evaluate(double *x, std::vector<double const *> const &args,
                    std::vector<std::vector<unsigned int> > const &dims) const;
                    
    void PenFactor(const double E, const double L, const double R,
   				const double mue,
				const double qQ,
				double& P, double& S) const;                     
                    
    void coul(int, double, double, double&, double&) const;
                
	bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) const;
	
	std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims,
    	
    	std::vector<double const *> const &values) const;
    	    	
	bool checkParameterValue(std::vector<double const *> const &args,
    	std::vector<std::vector<unsigned int> > const &dims) const;
}; 
    
}}//end namespaces

#endif /* SFACTOR3HEDP_H_ */

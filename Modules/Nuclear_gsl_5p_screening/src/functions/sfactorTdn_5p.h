#ifndef SFACTORTDN_5P_H_
#define SFACTORTDN_5P_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace nuclear {

/**
 * @short Astrophysical S Factors
 * sfactor3Hedp_5p returns the astrophysical S-factor as a function of energy.
 * <pre>
 * obsy1 = sfactorTdn_5p(obsx1, e1, gi, gf, ri, rf)
 * </pre>
 */

class sfactorTdn_5p : public ArrayFunction {
    public:
	sfactorTdn_5p();
	
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

#endif /* SFACTORTDN_5P_H_ */

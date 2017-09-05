#ifndef SFACTORTDNFAST_H_
#define SFACTORTDNFAST_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace bugs {

/**
 * @short Astrophysical S Factors
 * sfactorTdn returns the astrophysical S-factor as a function of energy.  It
 * returns a vector of S-factor predictions for a vector of energies.
 * <pre>
 * obsy1 = sfactorTdnF(obsx1, e1, gi, gf)
 * </pre>
 */

class sfactorTdn_fast : public ArrayFunction {
    public:
	sfactorTdn_fast();
	
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

#endif /* SFACTORTDNFAST_H_ */

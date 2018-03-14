#ifndef SFACTORTDN2_H_
#define SFACTORTDN2_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace nuclear {

/**
 * @short Astrophysical S Factors
 * sfactorTdn2 returns the astrophysical S-factor as a function of energy.
 * <pre>
 * obsy1 = sfactorTdn2(obsx1, e1, ex, gi, gf, ri, rf, ue)
 * </pre>
 */

class sfactorTdn2 : public ArrayFunction {
    public:
	sfactorTdn2();
	
    void evaluate(double *x, std::vector<double const *> const &args,
                    std::vector<std::vector<unsigned int> > const &dims) const;
                    
    void coul(int, double, double, double&, double&) const;
                
	bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) const;
	
	std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims,
    	
    	std::vector<double const *> const &values) const;
    	    	
	bool checkParameterValue(std::vector<double const *> const &args,
    	std::vector<std::vector<unsigned int> > const &dims) const;
}; 
    
}}//end namespaces

#endif /* SFACTORTDN2_H_ */
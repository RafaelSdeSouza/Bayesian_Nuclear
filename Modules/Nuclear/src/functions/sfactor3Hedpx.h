#ifndef SFACTOR3HEDPX_H_
#define SFACTOR3HEDPX_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace nuclear {

/**
 * @short Astrophysical S Factors; assumes B_c = Sc(E0), i.e., the
 * level shift is zero at the energy eigenvalue E0;
 *
 * sfactor3Hedpx returns the astrophysical S-factor as a function of energy.
 * <pre>
 * obsy1 = sfactor3Hedpx(obsx1, e0, ga, gb, ra, rb, ue)
 * </pre>
 */

class sfactor3Hedpx : public ArrayFunction {
    public:
	sfactor3Hedpx();
	
    void evaluate(double *x, std::vector<double const *> const &args,
             std::vector<std::vector<unsigned int> > const &dims) const;
                    
    void coul(int, double, double, double&, double&) const;
    
     void PenFactor(const double E, const double L, const double R,
   				const double mue, const double qQ, double& P, double& S) const;     
                
	bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) 
       const;
	
	std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > 
       const &dims,
    	
    	std::vector<double const *> const &values) const;
    	    	
	bool checkParameterValue(std::vector<double const *> const &args,
    	std::vector<std::vector<unsigned int> > const &dims) const;
}; 
    
}}//end namespaces

#endif /* SFACTOR3HEDPX_H_ */

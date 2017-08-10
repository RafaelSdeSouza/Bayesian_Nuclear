#ifndef SFACTORTDN_H_
#define SFACTORTDN_H_

#include <function/ArrayFunction.h>

namespace jags {
namespace nuclear {

/**
 * @short Astrophysical S Factors
 * sfactorTdn returns the astrophysical S-factor as a function of energy.  It
 * returns a 1000x2 table where column 1 is energy and column 2 is the S-factor.
 * It requires a fortran executable to perform the calculations.
 * <pre>
 * table = sfactorTdn(e1, gi, gf)
 * </pre>
 */

    class sfactorTdn : public ArrayFunction
    {
    public:
	sfactorTdn();
      void evaluate(double *x, std::vector<double const *> const &args,
                    std::vector<std::vector<unsigned int> > const &dims)
        const;
      bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) const;
      std::vector<unsigned int>
        dim(std::vector<std::vector<unsigned int> > const &dims,
            std::vector<double const *> const &values) const;
      bool checkParameterValue(std::vector<double const *> const &args,
                               std::vector<std::vector<unsigned int> > const &dims) const;
    }; 
}}

#endif /* SFACTORTDN_H_ */

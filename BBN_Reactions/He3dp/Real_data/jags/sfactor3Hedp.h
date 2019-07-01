
    extern "C" {
    #include "ArrayFunction.h"
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

};
    }
    
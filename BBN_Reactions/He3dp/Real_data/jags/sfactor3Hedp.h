
    extern "C" {
    #include </Users/Rafael/Documents/GitHub/UNC_Nuclear/JAGS-4.2.0/src/include/function/ArrayFunction.h>
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
    }
    
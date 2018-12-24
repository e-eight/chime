#include <complex>

namespace tdho
{
    class WaveFunction
    {
    public:
	WaveFunction();
	WaveFunction(double);
	WaveFunction(int, int, int, double);

	~WaveFunction();
	
	double norm_nl();
	double radial_nl(double);
	double theta_lm(double);
	std::complex<double> spherical_lm(double, double);

	double osc_b;
	int n, l, m;
    };
}

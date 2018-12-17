#include <iostream>
#include <cmath>

int delta(int a, int b)
{
	return a == b;
}

int main()
{
	int a = 1, b = 1;
	int d = delta(a, b);
	double s = std::sqrt(34);
	double ss = sqrt(34);
	std::cout << s << std::endl;
	std::cout << ss << std::endl;
}

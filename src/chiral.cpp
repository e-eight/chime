#include <iostream>
#include "chiral.h"

int main()
{
    chiral::Operator O(chiral::Operator::Names::identity);
    chiral::Operator::Names name = O.get_name();
    if (name == chiral::Operator::Names::identity)
	std::cout << "identity" << std::endl;
}

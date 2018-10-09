#ifndef __UTILITY_H__
#define __UTILITY_H__

inline int delta(int a, int b) { return a == b ? 1 : 0; };
double cg_coeff(int ja, int jb, int j, int mja, int mjb, int mj);

#define HEADING "ISU Two-body matrix elements\n"
#define OPERATOR "Operator: "
#define CONSTANTS "Relevant low energy constants: "
#define QNUMS "Quantum numbers: "

void print_preamble(char *operator, char *constants);

#endif

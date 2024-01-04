#include "thermofunctions.h"

double sroCorrectionFunction(double x, double a1, double b1, double a2){ 
	return a1 - a1*exp(b1/x) - a2/x + (a2/x)*exp(b1/x);
}

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <complex>
#include <math.h>

using namespace std;

class Eval {
  public:
    complex<double> eval(char *expr);
    complex<double> sum();
    complex<double> prod();
    complex<double> power();
    complex<double> basexp();
    complex<double> function();
    complex<double> number();
  private:
    void eatspace();
    char cur() { return *m_cptr; }
    void next();
    void expect(char ch);
    char *m_cptr;
};

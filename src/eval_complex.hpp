#ifndef EVAL_COMPLEX_HPP
#define EVAL_COMPLEX_HPP
#include <complex>

class Eval {
  public:
    std::complex<double> eval(char *expr);
    std::complex<double> sum();
    std::complex<double> prod();
    std::complex<double> power();
    std::complex<double> basexp();
    std::complex<double> function();
    std::complex<double> number();
  private:
    void eatspace();
    char cur() { return *m_cptr; }
    void next();
    void expect(char ch);
    char *m_cptr;
};
#endif /* EVAL_COMPLEX_HPP */

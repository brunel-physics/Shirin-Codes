#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <complex>
#include <math.h>
#include "eval_complex.hpp"
using namespace std;



void Eval::eatspace() {
  // überliest Spaces, Tabs usw.
  while(isspace(*m_cptr))
    m_cptr++;
}

void Eval::next() {
  // liest nächstes Zeichen aus Eingabe
  m_cptr++;
  eatspace();
}

void Eval::expect(char ch) {
  // erwartet vorgegebenes Zeichen aus Eingabe
  // Überspringen bei Übereinstimmung, sonst Exception
  if(cur() != ch)
    throw 0;
  next();
}

complex<double> Eval::eval(char *expr) {
  // wertet kompletten Ausdruck aus
  complex<double> res;
  m_cptr = expr;
  eatspace();
  res = sum();
  if(cur())
    throw 0;
  return res;
}

complex<double> Eval::sum() {
  // wertet Summe (evtl. mit nur 1 Summand) aus
  complex<double> res = prod();
  for(;;) {
    if(cur() == '+') {
      next();
      res += prod();
    }
    else if(cur() == '-') {
      next();
      res -= prod();
    }
    else
      break;
  }
  return res;
}

complex<double> Eval::prod() {
  // wertet Produkt (evtl. mit nur 1 Faktor) aus
  complex<double> res = power();
  for(;;) {
    if(cur() == '*') {
      next();
      res *= power();
    }
    else if(cur() == '/') {
      next();
      res /= power();
    }
    else
      break;
  }
  return res;
}

complex<double> Eval::power() {
  // wertet Potenz mit optionalem Vorzeichen (evtl. ohne Exponent) aus
  // ^ ist rechtsassoziativ und soll (wie in den meisten Programmiersprachen
  // mit Potenzoperator) stärker binden als das unäre + und -v (Vorzeichen),
  // d.h. -3^2 = -(3^2) = -9
  complex<double> res;
  switch(cur()) {
    case '+':
      next();
      res = power();
      break;
    case '-':
      next();
      res = -power();
      break;
    default:
      res = basexp();
      if(cur() == '^') {
        next();
        res = pow(res, power());
      }
  }
  return res;
}

complex<double> Eval::basexp() {
  // wertet Basis oder Exponent (Klammerausdruck oder Konstante) aus
  complex<double> res;
  if(cur() == '(') {
    next();
    res = sum();
    expect(')');
  }
  else {
    if(isalpha(cur()))
      res = function();
    else if(cur() == '.' || isdigit(cur()))
      res = number();
    else
      throw 0;
  }
  return res;
}

complex<double> complex_fabs(const complex<double>& num)
{
  double re = num.real();
  double im = num.imag();
  return sqrt(re*re + im*im);
}

complex<double> complex_ceil(const complex<double>& num)
{
  double re = ceil(num.real());
  double im = ceil(num.imag());
  return complex<double>(re, im);
}

complex<double> complex_floor(const complex<double>& num)
{
  double re = floor(num.real());
  double im = floor(num.imag());
  return complex<double>(re, im);
}

complex<double> re(const complex<double>& num)
{
  return num.real();
}

complex<double> im(const complex<double>& num)
{
  return num.imag();
}

complex<double> angle(const complex<double>& num)
{
  return atan2(num.imag(), num.real());
}

complex<double> conj(const complex<double>& num)
{
  return complex<double>(num.real(), -num.imag());
}

complex<double> complex_log2(const complex<double>& num)
{
  const complex<double> l2 = log(2.0);
  const complex<double> lnum = log(num);
  return lnum / l2;
}

double nthroot(double root, double num)
{
  return pow(num, 1/root);
}

complex<double> complex_round(const complex<double>& num)
{
  double re = round(num.real());
  double im = round(num.imag());
  return complex<double>(re, im);
}

complex<double> Eval::function() {
  // wertet symbolische Konstante oder Funktion mit einem oder zwei Argumenten
  // aus

  // Konstanten
  static struct {
    const char *name;
    complex<double> value;
  } consts[] = {
    { "pi",    3.14159265358979323846 },
    { "e",     2.71828182845904523536 },
    { "i",     1i },
    { "j",     1i }
  };

  // einstellige Funktionen
  static struct {
    const char *name;
    complex<double> (*func)(const complex<double>&);
  } funcs1[] = {
    { "sin",   sin },
    { "cos",   cos },
    { "tan",   tan },
    /*{ "asin",  asin },
    { "acos",  acos },
    { "atan",  atan },*/
    { "sinh",  sinh },
    { "cosh",  cosh },
    { "tanh",  tanh },
    /*{ "asinh", asinh },
    { "acosh", acosh },
    { "atanh", atanh },*/
    { "exp",   exp },
    { "log",   log },
    { "log2",  complex_log2 },
    { "log10", log10 },
    { "sqrt",  sqrt },
    { "fabs",  complex_fabs },
    { "ceil",  complex_ceil },
    { "floor", complex_floor },
    { "real", re },
    { "imag", im },
    { "angle", angle },
    { "conj", conj },
    { "round", complex_round }
  };

  // zweistellige Funktionen mit nur reellwertigen argumenten
  static struct {
    const char *name;
    double (*func)(double, double);
  } funcs2[] = {
    { "atan2", atan2 },
    { "hypot", hypot },
    { "fmod",  fmod },
    { "nthroot", nthroot }
  };

  char name[20];
  unsigned int i;

  i = 0;
  do {
    name[i++] = cur();
    next();
  } while(i < sizeof name - 1 && isalnum(cur()));
  if(i == sizeof name - 1)
    throw 1;
  name[i] = '\0';
  for(i = 0; i < sizeof consts / sizeof consts[0]; i++) {
    if(!strcmp(name, consts[i].name))
      return consts[i].value;
  }
  for(i = 0; i < sizeof funcs1 / sizeof funcs1[0]; i++) {
    if(!strcmp(name, funcs1[i].name)) {
      expect('(');
      complex<double> arg1 = sum();
      expect(')');
      return funcs1[i].func(arg1);
    }
  }
  for(i = 0; i < sizeof funcs2 / sizeof funcs2[0]; i++) {
    if(!strcmp(name, funcs2[i].name)) {
      expect('(');
      complex<double> arg1 = sum();
      expect(',');
      complex<double> arg2 = sum();
      expect(')');
      if(arg1.imag() != 0)
        throw 2;
      if(arg2.imag() != 0)
        throw 2;
      return funcs2[i].func(arg1.real(), arg2.real());
    }
  }

  throw 1;
}

complex<double> Eval::number() {
  // wertet Gleitkommazahl in Dezimaldarstellung aus
  double res = 0., p;
  int e, esign;

  while(isdigit(cur())) {
    // Vorkommastellen
    res = 10.0 * res + cur() - '0';
    next();
  }
  if(cur() == '.') {
    // Nachkommastellen
    next();
    p = 1.;
    while(isdigit(cur())) {
      p /= 10.;
      res += (double)(cur() - '0') * p;
      next();
    }
  }
  if(cur() == 'e' || cur() == 'E') {
    // Exponent
    next();
    esign = 1.;
    if(cur() == '+')
      next();
    else if(cur() == '-') {
      next();
      esign = -1.;
    }
    e = 0;
    while(isdigit(cur())) {
      e = 10 * e + cur() - '0';
      next();
    }
    res = res * pow(10., esign*e);
  }
  if(cur() == 'i' || cur() == 'j')
  {
    next();
    return complex<double>(0.0, res);
  }
  return complex<double>(res, 0.0);
}

//int main() {
//  static const char *message[] = {
//    "falsche Syntax",
//    "unbekannte Funktion oder Konstante",
//    "komplexes Argument"
//  };

  //Eval ev;
  //char line[1023];
  //complex<double> res;

 // fgets(line, sizeof line, stdin);
  //try {
    //res = ev.eval(line);
  //}
  //catch (int err) {
    //printf("Fehler: %s\n", message[err]);
    //return 1;
  //}
  //printf("%g   %g\n", res.real(), res.imag());
  //return 0;
//}

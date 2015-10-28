#if !defined(POLY_H)
/* ========================================================================
   grapt:
     Generates a LaTeX report of a user-defined graph property.

   poly.hpp:
     Basic polynomial class Poly, used for generating algebraic
     graphs.  
   
   Josh Tobin (tobinrj@tcd.ie), 2015
   ======================================================================== */

#define POLY_H

#include <iostream>
#include <vector>
#include <deque>
#include <utility>
#include <cstdlib>
#include <ctime>

// XXX get rid of this
using namespace std;

typedef vector<int> monomial;
typedef pair<int, monomial> term;

class poly
{
private:
  int n_vars;
  deque<term> terms;
public:

  poly(int _n_vars) : n_vars(_n_vars) {;}

  poly operator+(const term&);  
  poly operator+(const poly&);
  //poly operator-(const polt&);
  //poly operator*(const term&);

  poly& operator=(const poly& rhs);

  void print();
  int eval(const vector<int>&, int m); 
};

poly rand_poly(int d, int nv, int m);

#endif

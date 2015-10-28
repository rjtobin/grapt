/* ========================================================================
   grapt:
     Generates a LaTeX report of a user-defined graph property.

   poly.cpp:
     Basic polynomial class Poly, used for generating algebraic
     graphs.  
   
   Josh Tobin (tobinrj@tcd.ie), 2015
   ======================================================================== */

#include "poly.hpp"

static int mpow(int x, int p, int m)
{
  int ret = 1;
  for(int i=0; i<p; i++)
    ret = (ret*x)%m;
  return ret;
}

poly poly::operator+(const term& t)
{
  poly p = *this;
  
  monomial m = t.second;
  for(deque<term>::iterator i=p.terms.begin(); i!=p.terms.end(); i++)
  {
    if(i->second == m)
    {
      i->first += t.first;
      return p;
    }
  }

  p.terms.push_back(t);
  return p;
}

poly poly::operator+(const poly& p)
{
  poly ret = *this;

  // XXX oh god what is this (also we do this twice)
  deque<term> tmp_terms = p.terms;
  
  for(deque<term>::iterator i=tmp_terms.begin(); i!=tmp_terms.end(); i++)
  {
    ret = ret + *i;
  }

  return ret;
}

poly& poly::operator=(const poly& rhs)
{
  n_vars = rhs.n_vars;
  deque<term> tmp_terms = rhs.terms;
  
  terms.clear();
  for(deque<term>::iterator i=tmp_terms.begin(); i!=tmp_terms.end(); i++)
  {
    terms.push_back(*i);
  }

  return *this;
}

void poly::print()
{ 
  for(deque<term>::iterator i=terms.begin(); i!=terms.end(); i++)
  {
    monomial m = i->second;
    int coeff = i->first;
    if(coeff == 0)
      continue;
    if(coeff > 0)
    {
      cout << " + " << coeff << " ";
    }
    else
    {
      cout << " - " << -coeff << " ";
    }

    char c = 'w';

    for(int j=0; j<m.size(); j++)
    {
      c++;
      if(m[j] == 0)
        continue;
      cout << c << "^" << m[j] << ' ';
    }
  }
  cout << endl;
}

int poly::eval(const vector<int>& x, int m)
{
  int ret = 100*m;
  for(deque<term>::iterator i=terms.begin(); i!=terms.end(); i++)
  {
    int tmp = 1;
    monomial mon = i->second;
    for(int j=0; j<n_vars; j++)
    {
      tmp *= mpow(x[j],mon[j],m);
      tmp %= m;
    }
    tmp *= i->first;
    tmp = (tmp+100*m) % m;
    ret += tmp;
  }
  return (ret+100*m)%m;
}

/*int main()
{
  monomial m;
  m.resize(3);
  m[0] = m[1] = m[2] = 1;
  m[1] = 2;

  term t;
  t.first = 42;
  t.second = m;

  poly p(3);

  //p = p + t;
  //p.print();

  vector<int> x(3);
  x[0] = 1;
  x[1] = 2;
  x[2] = 17;

  p = rand_poly(1,3,17);

  p.print();
  
  cout << p.eval(x,1001) << endl;
  
  return 0;
  }*/

poly rand_poly(int d, int nv, int m)
{
  poly p(nv);

  monomial v(nv);

  for(int i=0; i<nv; i++)
    v[i] = 0;

  while(1)
  {
    bool done = false;
    for(int i=0; i<nv; i++)
    {
      if(v[i] < d)
      {
        v[i]++;
        done = true;
        break;
      }
      if(v[i] == d)
        v[i] = 0;
    }
    if(!done)
      break;
    int tmp = 0;
    for(int i=0; i<nv; i++)
      tmp += v[i];
    if(tmp > d)
      continue;

    term t;
    t.second = v;
    t.first = drand48() * m;
    p = p + t;
  }
  
  return p;
}

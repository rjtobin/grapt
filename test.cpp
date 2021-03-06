/* ========================================================================
   FILENAME $
   DESCRIPTION $
   PROJECT $
   Josh Tobin (tobinrj@tcd.ie), 2014
   ======================================================================== */

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "graph.hpp"

using namespace std;
using namespace arma;

#define EPS 0.0000001

double prec[500][500][500];

double extremal(int p, int q, int r)
{
  bool oob = false;
  if(p > 499 || q > 499 || r > 499)
    oob = true;
  if(!oob && prec[p][q][r] > -1.)
    return prec[p][q][r];

  mat A;
  int n = p+q+r;
  A.set_size(n,n);
  A.fill(1);
  for(int i=0; i<n; i++)
    A(i,i) = 0;
  for(int i=0; i<p; i++)
  {
    for(int j=n-1; j>=n-q; j--)
    {
      A(i,j) = A(j,i) = 0;
    }
  }

  vec eigval;
  mat eigvec;
  eig_sym(eigval,eigvec,A);

  double ans = eigval(n-1) + eigval(n-2);
  if(!oob)
    prec[p][q][r] = ans;  
  //prec[p][q][r] = eigval(n-1) + eigval(n-2);
  return ans;
}

bool w_test_graph(Graph* g, double* error)
{
  int n = g->getNumVertices();

  double t = g->spectrum(n-3);
  //double t2 = g->spectrum(2);
  //double t = (t1 > t2) ? t1 : t2;

  vector<int> deg(n);
  for(int i=0; i<n; i++)
    deg[i] = g->getDegree(i);

  sort(deg.begin(), deg.end());

  *error = deg[n-3] - t;
  if(*error < 0)
  {
    if(!isConnected(g))
    {
      *error = 1000;
      return true;
    }
    cout << deg[n-3] << ' ' << t << endl;
    return false;
  }
  return true;
    
}

bool _test_graph(Graph* g, double* error)
{
  int n = g->getNumVertices();

  double t1 = g->spectrum(n-1) + g->spectrum(n-2);

  double t2 = g->spectrum(0) + g->spectrum(n-2);

  double t = (t1 > t2) ? t1 : t2;

  //int p = 2. * n/7.;
  //int q = p;
  int r = 3. * n / 7.;
  int p = (n - r) * 0.5;
  int q = n - p - r;
  double bound = extremal(p,q,r);

  *error = bound + 0.000001 - t;

  //cout << bound << ' ' << t << endl;
  
  if(*error < 0)
    return false;
  return true;
}


bool test_graph(Graph* g, double* error)
{
  int n = g->getNumVertices();

  mat A1 = g->mAdjMatrix;
  mat A3 = g->mAdjMatrix;
  mat A5 = g->mAdjMatrix;
  mat A2 = g->mAdjMatrix;

  A3 = (A3*A3)*A3;

  A5 = (A3*A5)*A5;

  A2 = A1 * A1;
  
  double w1 = 0., w2 = 0., w3 = 0., w5 = 0.;

  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
    {
      w1 += A1(i,j);
      w2 += A2(i,j);
      w3 += A3(i,j);
      w5 += A5(i,j);
    }

  w1 /= (double) n;
  w3 /= (double) n;
  w5 /= (double) n;

  //*error =  2 * w5 * w1 - w3 * w3;
  *error =  n * w5 - w3 * w2;

  if(*error < -EPS)
    return false;
  return true;
}

bool ___test_graph(Graph* g, double* error)
{
  double total = 0;
  int n = g->getNumVertices();

  //if(g->spectrum(n-1) - g->spectrum(n-2) < EPS)
  //{
  //  *error = (double) n;
  //  return true;
  //}

  for(int i=n-1; i>=0; i--)
  {
    if(g->spectrum(i) > -EPS)
      total += pow(g->spectrum(i),2.);
  }

  if(total >= (double) (n-1) - EPS)
  {
    *error = total - (double) (n-1);
    return true;
  }

  if(!isConnected(g))
  {
    *error = ((double) n*n - (double) g->getNumEdges()) / ((double) (n*n) * 10.);
    return true;
  }

  *error = total - (double) (n-1);
  
  return false;
}

bool __test_graph(Graph* g, double* error)
{
  double E = 0;
  int n = g->getNumVertices();
  int m = g->getNumEdges();

  //if(g->spectrum(n-1) - g->spectrum(n-2) < EPS)
  //{
  //  *error = (double) n;
  //  return true;
  //}

  for(int i=n-1; i>=0; i--)
  {
    E += fabs(g->spectrum(i));
  }

  double rhs = 2. * (double) m / sqrt(g->spectrum(n-1));

  if(rhs - E > - EPS)
  {
    *error = rhs - E;
    return true;
  }

  if(!isConnected(g))
  {
    *error = ((double) n*n - (double) m) / ((double) (n*n) * 10.);
    return true;
  }

  *error = rhs - E;
  
  return false;
}

// elphick 1
bool el1_test_graph(Graph* g, double* error)
{
  double E = 0;
  int n = g->getNumVertices();
  int m = g->getNumEdges();

  //if(g->spectrum(n-1) - g->spectrum(n-2) < EPS)
  //{
  //  *error = (double) n;
  //  return true;
  //}

  for(int i=n-1; i>=0; i--)
  {
    E += fabs(g->spectrum(i));
  }

  double avg_second = 0.;
  for(int i=0; i<n; i++)
  {
    for(int j=i+1; j<n; j++)
    {
      if(g->isConnected(i,j))
      {
        avg_second += sqrt((double) (g->getDegree(i) * g->getDegree(j)));
      }
    }
  }
  
  double rhs = 2. * (double) (m*m) / avg_second;

  if(E - rhs > - EPS)
  {
    *error = E - rhs;
    return true;
  }

  if(!isConnected(g))
  {
    *error = 2. + ((double) n*n - (double) m) / ((double) (n*n));
    return true;
  }

  *error = E - rhs;
  
  return false;
}

bool randomGNPTest()
{
  Graph g(1);
  double error;
  for(int n=10; n<100; n+=5)
  {
    for(double p = 0.01; p<=1.; p+=0.1)
    {
      cout << "." << flush;
      for(int trial=0; trial<100; trial++)
      {
        //cout << trial << endl;
        randomGraphGNP(&g,p,n);
        if(!test_graph(&g,&error))
        {
          cout << "\nCounterexample: " << n << ' ' << p << endl;
          for(int i=0; i<n; i++)
          {
            //cout << g.spectrum(i) << ' ';
            cout << g.mAdjMatrix << endl;
          }
          //cout << endl;
          return false;
        }
      }
    }
  }
  cout << endl;
  return true;
}

static double norm(double* p, int n, double exp)
{
  double total = 0;
  for(int i=0; i<n; i++)
    total += pow(p[i],exp);
  return pow(total, 1. / exp);
}

bool randomGWTest()
{
  Graph g(1);
  double error;
  for(int n=30; n<100; n+=5)
  {
    cout << "trying n=" << n << endl;
    double p[n], p2[n];

    for(int i=0; i<n; i++)
      p[i] = drand48()*0.8 + 0.2;
    double norm_p = norm(p,n,2.);
    //for(int i=0; i<n; i++)
    //  p[i] /= norm_p;

    double min_v;
    int min_dir;
    int pm = 1;
    double inc = 0.1;

    for(int trials = 0; trials < 10000; trials++)
    {
      min_dir = -1;
      for(int plus_minus = -1; plus_minus <= 2; plus_minus += 2)
      {
        for(int i=0; i<n; i++)
        {
          for(int j=0; j<n; j++)
            p2[j] = p[j];
          p2[i] += inc * plus_minus;
          if(p2[i] <= 0 || p2[i] >= 1.)
            continue;
          //norm_p = norm(p2,n,2.);
          //for(int j=0; j<n; j++)
          //  p2[j] /= norm_p;
          double avg_error = 0;
          for(int t=0; t<200; t++)
          {
            randomGraphGW(&g,p2,n);
            if(!test_graph(&g,&error))
            {
             cerr << "Counterexample!" << endl;
             cerr << g.mAdjMatrix << endl;
             cerr << error << endl;
             return false;
            }
            avg_error += error;
          }
          error = avg_error / 200.;
          cout << i << ' ' << plus_minus << ' ' << error << endl;
          if(min_dir == -1 || error < min_v)
          {
            min_v = error;
            min_dir = i;
            pm = plus_minus;
          }
        }
      }
      p[min_dir] += inc * pm;
      if(min_v < 1.1)
        inc = 0.01;
      else if(min_v < 10.)
        inc = 0.05;
      else
        inc = 0.1;
      //norm_p = norm(p,n,2.);
      cout << "new p: (error = " << min_v << ")" << endl;
      for(int i=0; i<n; i++)
      {
        cout << p[i] << ' ';
        //p[i] /= norm_p;
      }
      cout << endl;
    }
  }
  return true;
}

bool randomGWJoinTest()
{
  Graph g(1), g1(1), g2(1);
  double error;
  for(int n=30; n<100; n+=10)
  {
    cout << "trying n=" << n << endl;
    double p[n], p2[n];

    for(int i=0; i<n; i++)
      p[i] = drand48()*0.8 + 0.2;
    double norm_p = norm(p,n,2.);
    //for(int i=0; i<n; i++)
    //  p[i] /= norm_p;

    double min_v;
    int min_dir;
    int pm = 1;
    double inc = 0.1;

    for(int trials = 0; trials < 10000; trials++)
    {
      min_dir = -1;
      for(int plus_minus = -1; plus_minus <= 2; plus_minus += 2)
      {
        for(int i=0; i<n; i++)
        {
          for(int j=0; j<n; j++)
            p2[j] = p[j];
          p2[i] += inc * plus_minus;
          if(p2[i] <= 0 || p2[i] >= 1.)
            continue;
          //norm_p = norm(p2,n,2.);
          //for(int j=0; j<n; j++)
          //  p2[j] /= norm_p;
          double avg_error = 0;
          for(int t=0; t<200; t++)
          {
            randomGraphGW(&g1,p2,n/2);
            randomGraphGW(&g2,p2+(n/2),n/2);

            joinGraphs(&g,&g1,&g2);
            
            if(!test_graph(&g,&error))
            {
             cerr << "Counterexample!" << endl;
             cerr << g.mAdjMatrix << endl;
             cerr << error << endl;
             return false;
            }
            avg_error += error;
          }
          error = avg_error / 200.;
          cout << i << ' ' << plus_minus << ' ' << error << endl;
          if(min_dir == -1 || error < min_v)
          {
            min_v = error;
            min_dir = i;
            pm = plus_minus;
          }
        }
      }
      p[min_dir] += inc * pm;
      if(min_v < 1.1)
        inc = 0.01;
      else if(min_v < 10.)
        inc = 0.05;
      else
        inc = 0.1;
      //norm_p = norm(p,n,2.);
      cout << "new p: (error = " << min_v << ")" << endl;
      for(int i=0; i<n; i++)
      {
        cout << p[i] << ' ';
        //p[i] /= norm_p;
      }
      cout << endl;
    }
  }
  return true;
}

bool kiteTest()
{
  Graph g(1);
  double error;
  for(int r=2; r<100; r++)
    for(int s=1; s<30; s++)
    {
      createKiteGraph(&g, r, s);
      if(!test_graph(&g, &error))
      {
        cout << "Counterexample!" << endl << endl;
        cout << g.mAdjMatrix << endl << endl;
        return false;
      }
      //cout << r << ' ' << s << ' ' << error << endl;
    }
  return true;
}

bool randomDiamTest()
{
  Graph g(1);
  double error;
  for(int n=50; n<1000; n+=50)
    for(int diam=10; diam<n/3; diam+=10)
    {
      for(double p=0.01; p<0.2; p+=0.01)
      {
        double av = 0.;
        for(int trial=0; trial<20; trial++)
        {
          randomGraphDiam(&g, diam, p, n);
          error = 0;
          if(!test_graph(&g, &error))
          {
            cout << "Counterexample!" << endl << endl;
            cout << g.mAdjMatrix << endl << endl;
            return false;
          }
          av += error;
        }
        av /= 20.;
        cout << n << ' ' << diam << ' ' << p << ' ' << av << endl;
      }
    }
  return true;
}

int main()
{
  for(int p=0; p<100; p++)
    for(int q=0; q<100; q++)
      for(int r=0; r<100; r++)
        prec[p][q][r] = -2.;
  srand48(10);

  cout << "Testing against kite graphs..." << endl << endl;

  if(!kiteTest())
    return 0;

  cout << "All passed." << endl << endl;

  cout << "Testing against high diameter..." << endl << endl;

  // if(!randomDiamTest())
  //  return 0;
  
  cout << "All passed." << endl << endl;
  
  cout << "Testing against GNP..." << endl << endl;
  
  if(!randomGNPTest())
    return 0;

  cout << "All passed." << endl << endl;

  cout << "Testing against dual GW..." << endl << endl;

  if(!randomGWJoinTest())
    return 0;

  cout << "All passed." << endl << endl;
  
  cout << "Testing against GW..." << endl << endl;

  if(!randomGWTest())
    return 0;

  cout << "All passed." << endl;
  
  return 0;
}

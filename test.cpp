/* ========================================================================
   FILENAME $
   DESCRIPTION $
   PROJECT $
   Josh Tobin (tobinrj@tcd.ie), 2014
   ======================================================================== */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "graph.hpp"

using namespace std;

#define EPS 0.0000001

bool test_graph(Graph* g, double* error)
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
    *error = (double) n*n - (double) g->getNumEdges();
    return true;
  }

  *error = total - (double) (n-1);
  
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
      cout << ".";
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
        }
      }
    }
  }
  cout << endl;
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
  for(int n=10; n<100; n+=5)
  {
    cout << "trying n=" << n << endl;
    double p[n], p2[n];

    for(int i=0; i<n; i++)
      p[i] = drand48()*0.5 + 0.5;
    double norm_p = norm(p,n,2.);
    //for(int i=0; i<n; i++)
    //  p[i] /= norm_p;

    double min_v;
    int min_dir;
    int pm = 1;

    for(int trials = 0; trials < 1000; trials++)
    {
      min_dir = -1;
      for(int plus_minus = -1; plus_minus <= 2; plus_minus += 2)
      {
        for(int i=0; i<n; i++)
        {
          for(int j=0; j<n; j++)
            p2[j] = p[j];
          p2[i] += 0.01 * plus_minus;
          if(p2[i] <= 0 || p2[i] >= 1.)
            continue;
          //norm_p = norm(p2,n,2.);
          //for(int j=0; j<n; j++)
          //  p2[j] /= norm_p;
          double avg_error = 0;
          for(int t=0; t<1000; t++)
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
          error = avg_error / 1000.;
          cout << i << ' ' << plus_minus << ' ' << error << endl;
          if(min_dir == -1 || error < min_v)
          {
            min_v = error;
            min_dir = i;
            pm = plus_minus;
          }
        }
      }
      p[min_dir] += 0.01 * pm;
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

int main()
{
  srand48(10);

  cout << "Testing against GNP..." << endl << endl;
  
  //randomGNPTest();

  cout << "Done!" << endl << endl;

  cout << "Testing against GW..." << endl << endl;

  randomGWTest();

  cout << "Done!" << endl;
  
  return 0;
}

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

double statistic(Graph* g)
{
  int n = g->getNumVertices();
  int m = g->getNumEdges();
  double rho = g->spectrum(n-1);
  
  return rho - 2. * m / n;
}

void createPineapple(Graph* g, int n, int q)
{
  g->setNumVertices(n);
  for(int i=0; i<q; i++)
    for(int j=i+1; j<q; j++)
      g->addEdge(i,j);
  for(int i=q; i<n; i++)
    g->addEdge(q-1,i);
}

bool test_graph(Graph* g, double* error, double bound)
{
  int n = g->getNumVertices();
  int m = g->getNumEdges();
  
  double stat = statistic(g);

  if(!isConnected(g))
  {
    *error = 2. + ((double) n*n - (double) m) / ((double) (n*n));
    return true;
  }

  
  if(bound - stat > - EPS)
  {
    *error = bound - stat;
    return true;
  }

  *error = bound - stat;
  
  return false;
}

bool randomGNPTest()
{
  Graph g(1), rG(1);
  double error;
  for(int n=10; n<100; n+=5)
  {
    createPineapple(&rG,n,(int)(n/2) + 1);
    double bound = statistic(&rG);
    for(double p = 0.01; p<=1.; p+=0.1)
    {
      cout << ".";
      for(int trial=0; trial<100; trial++)
      {
        //cout << trial << endl;
        randomGraphGNP(&g,p,n);
        if(!test_graph(&g,&error,bound))
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
  Graph g(1), rG(1);
  double error;
  for(int n=15; n<100; n+=5)
  {
    cout << "trying n=" << n << endl;
    double p[n], p2[n];

    createPineapple(&rG, n, (int)(n/2) + 1);
    double bound = statistic(&rG);
     
    for(int i=0; i<n; i++)
      p[i] = drand48()*0.8 + 0.2;
    double norm_p = norm(p,n,2.);
    //for(int i=0; i<n; i++)
    //  p[i] /= norm_p;

    double min_v;
    int min_dir;
    int pm = 1;
    double inc = 0.1;

    for(int trials = 0; trials < 1000; trials++)
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
            if(!test_graph(&g,&error,bound))
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
  Graph g(1), rG(1);
  double error, comparison;
  for(int r=2; r<100; r++)
    for(int s=1; s<30; s++)
    {
      if(r+s-1 < 10)
        continue;
      createKiteGraph(&g, r, s);
      createPineapple(&rG, r+s-1, ((int) ((r+s-1)/2)) + 1);
      comparison = statistic(&rG);
      if(!test_graph(&g, &error, comparison))
      {
        cout << "Counterexample!" << endl << endl;
        cout << g.mAdjMatrix << endl << endl;
        return false;
      }
    }
  return true;
}

int main()
{
  srand48(10);

  Graph pine(1), kite(1);
  createPineapple(&pine,10,6);
  createKiteGraph(&kite,3,8);

  cout << statistic(&pine) << ' ' << statistic(&kite) << endl;
  
  cout << "Testing against kite graphs..." << endl << endl;

  if(!kiteTest())
    return 0;

  cout << "All passed." << endl << endl;
  
  cout << "Testing against GNP..." << endl << endl;
  
  if(!randomGNPTest())
    return 0;

  cout << "All passed." << endl << endl;

  cout << "Testing against GW..." << endl << endl;

  if(!randomGWTest())
    return 0;

  cout << "All passed." << endl;
  
  return 0;
}

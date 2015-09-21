/* ========================================================================
   FILENAME $
   DESCRIPTION $
   PROJECT $
   Josh Tobin (tobinrj@tcd.ie), 2014
   ======================================================================== */

#include <iostream>

#include "graph.hpp"

using namespace std;
using namespace arma;

int diam_main()
{
  Graph g(10);

  for(int diam=20; diam<50; diam+=10)
  {
    randomGraphDiam(&g, diam, 0.5, 1000);
    cout << "rg" << endl << endl << flush;
    cout << g.getDiam() << endl << flush;
  }
  
  return 0;
}

int _main()
{
  Graph g(5);
  int i,j;
  do
  {
    cin >> i >> j;
    g.addEdge(i,j);
  } while(i != -1);

  cout << g.mAdjMatrix << endl << endl;
  cout << g.getLaplacian() << endl << endl;
  cout << g.getNormLaplacian() << endl << endl;
  
  return 0;
}

int main()
{
  srand48(time(NULL));
  
  Graph g(1);

  randomGraphPoly(&g, 29, 2, 2, 3);

  cout << g.getNumEdges() << endl;
}


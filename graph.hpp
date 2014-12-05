#if !defined(GRAPH_H)
/* ========================================================================
   FILENAME $
   DESCRIPTION $
   PROJECT $
   Josh Tobin (tobinrj@tcd.ie), 2014
   ======================================================================== */

#define GRAPH_H

#include <armadillo>
#include <cstdlib>
#include <deque>

class Graph
{
public:
  Graph(int n);
  Graph(arma::mat* adjacency);
  ~Graph();
  
  void setNumVertices(int n);
   
  bool isConnected(int i, int j);
  void addEdge(int i, int j);

  double spectrum(int i);
  double eigenvector(int i, int j);

  int getNumVertices();
  int getNumEdges();

  int getDegree(int i);
  
//private:
  arma::mat mAdjMatrix;
  int mN;
  int mE;
  int* mDeg;
  double* mSpectrum;
  arma::mat mEigenvectors;
  
  bool mSpectOutdated;
  void mGenerateSpectrum();
};

void randomGraphGNP(Graph* g, double p, int n);
void randomGraphGW(Graph* g, double* p, int n);

void createKiteGraph(Graph* g, int pathSize, int cliqueSize);

bool isConnected(Graph* g);

#endif


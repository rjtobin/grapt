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

  int getDiam();
  
//private:
  arma::mat mAdjMatrix;
  int mN;
  int mE;
  int* mDeg;
  double* mSpectrum;
  int mDiam;
  arma::mat mEigenvectors;
  
  bool mSpectOutdated;
  bool mDiamOutdated;
  
  void mGenerateSpectrum();
};

void randomGraphGNP(Graph* g, double p, int n);
void randomGraphGW(Graph* g, double* p, int n);
void randomGraphDiam(Graph* g, int diam, double p, int n);

void joinGraphs(Graph* r, Graph* g1, Graph* g2);

void createKiteGraph(Graph* g, int pathSize, int cliqueSize);

bool isConnected(Graph* g);

#endif


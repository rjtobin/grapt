#if !defined(GRAPH_H)
/* ========================================================================
   grapt:
     Generates a LaTeX report of a user-defined graph property.

   graph.hpp:  Defines the Graph class 
   
   Josh Tobin (tobinrj@tcd.ie), 2015
   ======================================================================== */

#define GRAPH_H

#include <armadillo>
#include <cstdlib>
#include <deque>

#include "poly.hpp"

class Graph
{
public:
  Graph(int n=1);
  Graph(arma::mat* adjacency);
  Graph(const Graph& copy_from);
  Graph(char* g6_graph);  
  ~Graph();

  void from_g6(const char* d);
  void to_g6(char* d);
  
  Graph& operator=(const Graph& rhs);
  bool operator<(const Graph& rhs) const;
  bool operator>(const Graph& rhs) const;
  bool operator==(const Graph& rhs) const;
  
  void setNumVertices(int n);
   
  bool isConnected(int i, int j);
  void addEdge(int i, int j);
  void deleteEdge(int i, int j);

  double spectrum(int i);
  double eigenvector(int i, int j);
  double mainAngle(int i);

  int getNumVertices();
  int getNumEdges();

  int getDegree(int i);

  int getDiam();

  arma::mat& getNormLaplacian();
  arma::mat& getLaplacian();
  
//private:
  arma::mat mAdjMatrix;

private:  
  arma::mat mNormLaplacian;
  arma::mat mLaplacian;
  
  int mN;
  int mE;
  int* mDeg;
  double* mSpectrum;
  int mDiam;
  arma::mat mEigenvectors;
  
  bool mSpectOutdated;
  bool mDiamOutdated;
  bool mNLOutdated;
  
  void mGenerateSpectrum();
  void mGenerateNL();
};

void randomGraphGNP(Graph* g, double p, int n);
void randomGraphGW(Graph* g, double* p, int n);
void randomGraphDiam(Graph* g, int diam, double p, int n);
void randomGraphPoly(Graph* g, int q, int t, int r, int d);

void joinGraphs(Graph* r, Graph* g1, Graph* g2);

void createKiteGraph(Graph* g, int pathSize, int cliqueSize);

bool isConnected(Graph* g);

#endif


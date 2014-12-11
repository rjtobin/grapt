/* ========================================================================
   FILENAME $
   DESCRIPTION $
   PROJECT $
   Josh Tobin (tobinrj@tcd.ie), 2014
   ======================================================================== */
#include "graph.hpp"

#define EPS 0.0000001

using namespace arma;
using namespace std;

static inline int min(int x, int y)
{
  return (x<y) ? x : y;
}

Graph::Graph(int n)
{
  mN = n;
  mE = 0;
  mAdjMatrix.set_size(n,n);
  mAdjMatrix.fill(0);
  mSpectrum = 0;
  mSpectOutdated = true;
  mDiamOutdated = true;
  mDeg = 0;
}

Graph::Graph(mat* adjacency)
{
  mAdjMatrix = *adjacency;
  mN = mAdjMatrix.n_cols;
  mSpectrum = 0;
  mSpectOutdated = true;
  mDiamOutdated = true;

  mDeg = new int[mN];
  
  mE = 0;
  for(int i=0; i<mN; i++)
  {
    mDeg[i] = 0;
    for(int j=i; j<mN; j++)
    {
      if(mAdjMatrix(i,j) > 0.5)
      {
        mDeg[i]++;
        mE++;
      }
    }
  }
}

Graph::~Graph()
{
  if(mSpectrum)
    delete[] mSpectrum;
  if(mDeg)
    delete[] mDeg;
}

void Graph::setNumVertices(int n)
{
  mN = n;
  mE = 0;
  mAdjMatrix.resize(n,n);
  mAdjMatrix.fill(0);
  mSpectOutdated = true;
  mDiamOutdated = true;
  if(mDeg)
    delete[] mDeg;
  mDeg = new int[n];
  for(int i=0; i<n; i++)
    mDeg[i] = 0;
}

bool Graph::isConnected(int i, int j)
{
  if(i < 0 || j < 0 || i >= mN || j >= mN)
    return false;
  return (mAdjMatrix(i,j) > EPS);
}

void Graph::addEdge(int i, int j)
{
  if(i < 0 || j < 0 || i >= mN || j >= mN)
    return;
  if(mAdjMatrix(i,j) < EPS)
  {
    mAdjMatrix(i,j) = mAdjMatrix(j,i) = 1.;
    mE++;
    mDeg[i]++;
    mDeg[j]++;
    mSpectOutdated = true;
    mDiamOutdated = true;
  }
}

double Graph::spectrum(int i)
{
  if(i<0 || i>=mN)
    return 0;
  if(mSpectOutdated)
    mGenerateSpectrum();
  return mSpectrum[i];
}

double Graph::eigenvector(int i, int j)
{
  if(i<0 || j<0 || i>=mN || j>=mN) 
    return 0;
  if(mSpectOutdated)
    mGenerateSpectrum();
  return mEigenvectors(i,j);
}

int Graph::getNumVertices()
{
  return mN;
}

int Graph::getNumEdges()
{
  return mE;
}

int Graph::getDegree(int i)
{
  return mDeg[i];
}

int Graph::getDiam()
{
  if(!mDiamOutdated)
    return mDiam;

  mat D, Dn;
  //D.set_size(mN,mN);
  D = mAdjMatrix;
  for(int i=0; i<mN; i++)
    for(int j=0; j<mN; j++)
      if(i != j && D(i,j) < 0.5)
        D(i,j) = mN*mN*mN; // infinity XXX: improve this 

  Dn = D;
  for(int k=0; k<mN; k++)
  {
    for(int i=0; i<mN; i++)
    {
      for(int j=0; j<mN; j++)
      {
        Dn(i,j) = min(D(i,j), D(i,k) + D(k,j));
      }
    }
    D = Dn;
  }
  mDiam = 0;
  for(int i=0; i<mN; i++)
    for(int j=0; j<mN; j++)
      if(D(i,j) > mDiam)
        mDiam = D(i,j);
  mDiamOutdated = false;
  return mDiam;
}

void Graph::mGenerateSpectrum()
{
  vec eigval;
  eig_sym(eigval, mEigenvectors, mAdjMatrix);

  if(mSpectrum)
    delete[] mSpectrum;
  mSpectrum = new double[mN];
  
  for(int i=0; i<mN; i++)
    mSpectrum[i] = eigval(i);
  mSpectOutdated = false;
}

void randomGraphGNP(Graph* g, double p, int n)
{
  g->setNumVertices(n);
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
    {
      if(drand48() < p)
        g->addEdge(i,j);
    }
}

void randomGraphGW(Graph* g, double* p, int n)
{
  g->setNumVertices(n);
  
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
    {
      if(drand48() < p[i]*p[j])
        g->addEdge(i,j);
    }
}

void createKiteGraph(Graph* g, int r, int s)
{
  int n = r+s-1;
  g->setNumVertices(r+s-1);

  int cliqueStart = r;
  for(int i=cliqueStart; i<n; i++)
    for(int j=i+1; j<n; j++)
      g->addEdge(i,j);

  for(int i=0; i<cliqueStart; i++)
    g->addEdge(i,i+1);
    
}

bool isConnected(Graph* g)
{
  int n = g->getNumVertices();

  bool visited[n];
  int component[n];
  for(int i=0; i<n; i++)
  {
    visited[i] = false;
    component[i] = i;
  }
  
  deque<int> next;
  next.push_back(0);
  visited[0] = true;
  
  while(!next.empty())
  {
    int cur = next.back();
    next.pop_back();

    int min_v = n+1, min_n;
    for(int i=0; i<n; i++)
      if(g->isConnected(i,cur))
      {
        if(component[i] < min_v)
        {
          min_n = i;
          min_v = component[i];
        }
        if(!visited[i])
        {
          next.push_back(i);
          visited[i] = true;
        }
      }
    if(min_v < n+1)
      component[cur] = min_v;
  }

  for(int i=0; i<n; i++)
    if(visited[i] != true)
      return false;
  return true;
}

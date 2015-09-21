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
  mLaplacian.set_size(n,n);
  mLaplacian.fill(0);
  mSpectrum = 0;
  mSpectOutdated = true;
  mDiamOutdated = true;
  mNLOutdated = true;

  mDeg = new int[n];
}

Graph::Graph(const Graph& copy_from)
{
  mN = copy_from.mN;
  mE = copy_from.mE;
  mAdjMatrix = copy_from.mAdjMatrix;
  mLaplacian = copy_from.mLaplacian;
  mSpectrum = 0;
  mSpectOutdated = true;
  mDiamOutdated = true;
  mNLOutdated = true;

  mDeg = new int[mN];
  for(int i=0; i<mN; i++)
    mDeg[i] = copy_from.mDeg[i];
}

Graph& Graph::operator=(const Graph& rhs)
{
  mN = rhs.mN;
  mE = rhs.mE;
  mAdjMatrix = rhs.mAdjMatrix;
  mLaplacian = rhs.mLaplacian;
  if(mSpectrum)
    delete[] mSpectrum;
  mSpectrum = 0;
  mSpectOutdated = true;
  mDiamOutdated = true;
  mNLOutdated = true;

  if(mDeg)
    delete[] mDeg;
  mDeg = new int[mN];
  for(int i=0; i<mN; i++)
    mDeg[i] = rhs.mDeg[i];
}


Graph::Graph(mat* adjacency)
{
  mAdjMatrix = *adjacency;
  mN = mAdjMatrix.n_cols;
  mSpectrum = 0;
  mSpectOutdated = true;
  mDiamOutdated = true;
  mNLOutdated = true;

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

  mLaplacian.fill(0);
  for(int i=0; i<mN; i++)
  {
    mLaplacian(i,i) = mDeg[i];
  }
  mLaplacian = mLaplacian - mAdjMatrix;
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

  mLaplacian.resize(n,n);
  mLaplacian.fill(0);
  
  mSpectOutdated = true;
  mDiamOutdated = true;
  mNLOutdated = true;
  
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
    mLaplacian(i,j) = mLaplacian(j,i) = -1.;
    mLaplacian(i,i) += 1.;
    mLaplacian(j,j) += 1.;
    mE++;
    mDeg[i]++;
    mDeg[j]++;
    mSpectOutdated = true;
    mDiamOutdated = true;
    mNLOutdated = true;
  }
}

void Graph::deleteEdge(int i, int j)
{
  if(i < 0 || j < 0 || i >= mN || j >= mN)
    return;
  if(mAdjMatrix(i,j) > 0.5)
  {
    mAdjMatrix(i,j) = mAdjMatrix(j,i) = 0.;
    mLaplacian(i,j) = mLaplacian(j,i) = 0.;
    mLaplacian(i,i) -= 1.;
    mLaplacian(j,j) -= 1.;
    mE--;
    mDeg[i]--;
    mDeg[j]--;
    mSpectOutdated = true;
    mDiamOutdated = true;
    mNLOutdated = true;
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

  //std::cerr << mN << ' ' << endl;
 
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

arma::mat& Graph::getNormLaplacian()
{
  if(mNLOutdated)
    mGenerateNL();
  return mNormLaplacian;
}

arma::mat& Graph::getLaplacian()
{
  return mLaplacian;
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

void Graph::mGenerateNL()
{
  mNormLaplacian.set_size(mN,mN);
  mNormLaplacian.fill(0);
  for(int i=0; i<mN; i++)
  {
    for(int j=0; j<mN; j++)
    {
      if(i==j)
        mNormLaplacian(i,j) = 1.;
      else if(mAdjMatrix(i,j) > 0.5)
        mNormLaplacian(i,j) = 1. / sqrt(mDeg[i] * mDeg[j]);
    }
  }
  mNLOutdated = false;
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

void randomGraphDiam(Graph* g, int diam, double p, int n)
{
  g->setNumVertices(n);
  
  int vert[diam];
  for(int i=0; i<diam; i++)
    vert[i] = 1;

  int left = n-diam;
  while(left--)
  {
    int index = drand48() * diam;
    vert[index]++;
  }

  g->setNumVertices(n);

  int current_vertex = 0;
  
  for(int i=0; i<diam; i++)
  {
    for(int j=0; j<vert[i]; j++)
      for(int k=j+1; k<vert[i]; k++)
        if(drand48() < p)
          g->addEdge(current_vertex+j, current_vertex+k);

    if(i<diam-1)
    {
      int v1 = vert[i], v2 = vert[i+1];
      
      for(int j=0; j<v1; j++)
        for(int k=0; k<v2; k++)
          if(drand48() < p)
            g->addEdge(current_vertex + j, current_vertex + v1 + k);
    }
    
    current_vertex += vert[i];
  }  
}

void randomGraphPoly(Graph* g, int q, int t, int r, int d)
{
  int n = 2;
  for(int i=0; i<t; i++)
    n*=q;
  g->setNumVertices(n);

  poly* p[r];
  for(int i=0; i<r; i++)
  {
    p[i] = new poly(2*t);
    *(p[i]) = rand_poly(d,2*t,q);
    p[i]->print();
  }

  vector<int> v(2*t);
  for(int i=0; i<2*t; i++)
    v[i] = 0;

  while(1)
  {
    bool done = true;

    for(int i=0; i<2*t; i++)
    {
      if(v[i] < q-1)
      {
        v[i]++;
        done = false;
        break;
      }
      else
        v[i] = 0;
    }

    bool connect = true;
    for(int j=0; j<r; j++)
    {
      if( (p[j])->eval(v,q) != 0)
      {
        connect = false;
        break;
      }
    }

    if(connect)
    {
      int index1 = 0, index2 = 0;
      int qpow = 1;
      for(int i=0; i<t; i++)
      {
        index1 += v[i] * qpow;
        index2 += v[t+i] * qpow;
        qpow *= q;
      }
      g->addEdge(index1,index2);
      //for(int i=0; i<2*t; i++)
      //  cout << v[i] << ' ';
      //cout << endl;
    }
    
    if(done)
      break;
  }
  
  for(int i=0; i<r; i++)
    delete (p[i]);
}

void joinGraphs(Graph* r, Graph* g1, Graph* g2)
{
  int n = g1->getNumVertices();
  int m = g2->getNumVertices();

  r->setNumVertices(n+m);

  for(int i=0; i<n; i++)
  {
    for(int j=i+1; j<n; j++)
    {
      if(g1->isConnected(i,j))
        r->addEdge(i,j);
    }
  }

  for(int i=0; i<m; i++)
  {
    for(int j=0; j<m; j++)
    {
      if(g2->isConnected(i,j))
        r->addEdge(n+i,n+j);
    }
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


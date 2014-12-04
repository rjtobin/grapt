/* ========================================================================
   FILENAME $
   DESCRIPTION $
   PROJECT $
   Josh Tobin (tobinrj@tcd.ie), 2014
   ======================================================================== */

#include <iostream>
#include <cmath>
#include "graph.hpp"

using namespace std;

int main()
{
  int k = 101;
  int m = 30;
  int start_of_clique = ((k-1)/2) * m;
  
  Graph G(start_of_clique + k + 1);
  
  for(int i=0; i<m; i++)
  {
    //int prev = (i-1 + m) % m;
    int next = (i+1 + m) % m;

    for(int j=0; j<(k-1)/2; j++)
      for(int l=0; l<(k-1)/2; l++)
      {
        //G.addEdge(prev * ((k-1)/2) + j, i * ((k-1)/2) + l);
        if(!(i == 0 && j == 0 && l == 0)) 
          G.addEdge(next * ((k-1)/2) + j, i * ((k-1)/2) + l);
      }
  }

  G.addEdge(start_of_clique, 0);
  G.addEdge(start_of_clique + 1, (k-1)/2);
  
  for(int i=0; i<k+1; i++)
    for(int j=i+1; j<k+1; j++)
      if(i != 0 || j != 1)
        G.addEdge(start_of_clique + i, start_of_clique + j);

  double total = 0.;
  for(int i=0; i<start_of_clique+k; i++)
  {
    total += pow(G.eigenvector(i,start_of_clique+k),1);
    //cout << i << ' ' << 
  }

  total *= total;
  total /= (double) (start_of_clique+k+1); 
  
  cout << G.spectrum(start_of_clique+k) << ' ' << total << endl; 

   
  //cout << G.mAdjMatrix << endl;
  
  return 0;
}

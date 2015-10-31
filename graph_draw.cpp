/* ========================================================================
   grapt:
     Generates a LaTeX report of a user-defined graph property.

   graph_draw.cpp:  Force-directed graph drawing functions. 
   
   Josh Tobin (tobinrj@tcd.ie), 2015
   ======================================================================== */

#include "graph_draw.hpp"

#define EPS 1e-9

void draw_graph(CDrawing& drawing, Graph& g, double x_size, double y_size)
{
  int n = g.getNumVertices();

  CPoint* pts[n];
  for(int i=0; i<n; i++)
  {
    pts[i] = new CPoint(drand48() * (double) x_size, drand48() * (double) y_size);
    drawing.addShape(pts[i]);
  }

  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
      if(g.isConnected(i,j))
        drawing.addShape(new CLine(*pts[i],*pts[j]));
}

// Force algorithm of Eades (1984)
void force_draw(CDrawing& drawing, Graph& g, double x_size, double y_size)
{
  int n = g.getNumVertices();

  double c1 = 2., c2 = 1., c3 = 1., c4 = .1;
  int n_iter = 400;
  
  CPoint* pts[n];
  CPoint* forces[n];
  for(int i=0; i<n; i++)
  {
    forces[i] = new CPoint(0,0);
    pts[i] = new CPoint(drand48() * (double) (x_size / 1.2), drand48() * (double) (y_size / 1.2));
  }

  if(n > 40)
    n_iter = 0;

  for(int i=0; i<n_iter; i++)
  {
    for(int j=0; j<n; j++)
    {
      double fx = 0.,  fy = 0., dist, fc;
      for(int k=0; k<n; k++)
      {
        if(j==k)
          continue;
        dist = sqrt(pow(pts[j]->x - pts[k]->x,2.) + pow(pts[j]->y - pts[k]->y,2.));
        if(g.isConnected(j,k))
          fc = - c1 * log(dist / c2);
        else
          fc = c3 / pow(dist,2.);
        
        fc *= c4;        

        fx += fc * (pts[j]->x - pts[k]->x) / dist;
        fy += fc * (pts[j]->y - pts[k]->y) / dist;
      }

      pts[j]->x += fx;
      pts[j]->y += fy;

    }
  }

  double max_x = pts[0]->x, min_x = pts[0]->x;
  double max_y = pts[0]->x, min_y = pts[0]->x;
  
  for(int i=1; i<n; i++)
  {
    if(pts[i]->x > max_x)
      max_x = pts[i]->x;
    if(pts[i]->y > max_y)
      max_y = pts[i]->y;
    if(pts[i]->x < min_x)
      min_x = pts[i]->x;
    if(pts[i]->y < min_y)
      min_y = pts[i]->y;    
  }

  double width = max_x - min_x;
  double height = max_y - min_y;

  if(width < EPS)
    width = 1.;
  if(height < EPS)
    height = 1.;

  
  // If the graph is already within the specified size,
  // then trust the dimensions resulting from the algorithm
  if(width < x_size)
    width = x_size = 1.;
  if(height < y_size)
    height = y_size = 1.;
  
  
  for(int i=0; i<n; i++)
  {
     pts[i]->x = ((pts[i]->x - min_x) / width) * x_size;
     pts[i]->y = ((pts[i]->y - min_y) / height) * y_size;    
  }
  
  for(int i=0; i<n; i++)
  {
    drawing.addShape(pts[i]);
  }
  
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
      if(g.isConnected(i,j))
        drawing.addShape(new CLine(*pts[i],*pts[j]));  
  
  for(int i=0; i<n; i++)
    delete forces[i];
}


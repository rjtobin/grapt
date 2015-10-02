/* ========================================================================
   FILENAME $
   DESCRIPTION $
   PROJECT $
   Josh Tobin (tobinrj@tcd.ie), 2014
   ======================================================================== */

#include "graph_draw.hpp"

void draw_graph(CDrawing& drawing, Graph& g, int x_size, int y_size)
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
void force_draw(CDrawing& drawing, Graph& g, int x_size, int y_size)
{
  int n = g.getNumVertices();

  double c1 = 2., c2 = 1., c3 = 1., c4 = .1;
  int n_iter = 400;

  /*if(n > 30)
  {
    c3 /= 2.;
    c4 / 1.5;
    }*/
  
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

//      if(pts[j]->x < -5. || pts[j]->x > x_size || pts[j]->y < 5.
//         || pts[j]->y > y_size)
//      {
//        c3 = c3*0.99;
//        n_iter = 1;
//      }
    }
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


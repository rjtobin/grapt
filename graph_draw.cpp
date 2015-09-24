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

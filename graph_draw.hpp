#if !defined(GRAPH_DRAW_H)
/* ========================================================================
   grapt:
     Generates a LaTeX report of a user-defined graph property.

   graph_draw.hpp:  Force-directed graph drawing functions. 
   
   Josh Tobin (tobinrj@tcd.ie), 2015
   ======================================================================== */

#define GRAPH_DRAW_H

// XXX get rid of this!
#include <cstdlib>

#include "../clatex/clatex.hpp"
#include "graph.hpp"

void draw_graph(CDrawing& drawing, Graph& g, int x_size, int y_size);

void force_draw(CDrawing& drawing, Graph& g, int x_size, int y_size);

#endif

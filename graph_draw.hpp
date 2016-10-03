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

#include <emmintrin.h>
#include <pmmintrin.h>

#include "../clatex/clatex.hpp"
#include "graph.hpp"

void draw_graph(CDrawing& drawing, Graph& g, double x_size, double y_size);

void force_draw(CDrawing& drawing, Graph& g, double x_size, double y_size, int num_iterations, std::string* labels=0);
void sse_force_draw(CDrawing& drawing, Graph& g, double x_size, double y_size, int num_iterations);

#endif

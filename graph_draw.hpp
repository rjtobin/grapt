#if !defined(GRAPH_DRAW_H)
/* ========================================================================
   FILENAME $
   DESCRIPTION $
   PROJECT $
   Josh Tobin (tobinrj@tcd.ie), 2014
   ======================================================================== */

#define GRAPH_DRAW_H

// XXX get rid of this!
#include <cstdlib>

#include "../clatex/clatex.hpp"
#include "graph.hpp"

void draw_graph(CDrawing& drawing, Graph& g, int x_size, int y_size);

#endif

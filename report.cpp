/* ========================================================================
   grapt:
     Generates a LaTeX report of a user-defined graph property.

   report.cpp:
     Generates a report of the property defined by the function
     TestProperty(Graph* g)

   NOTE:
     To use this, the directory containing the executable
     must contain directories trees/, connected/ and reg/,
     which contain .g6 datafiles with all of these graph.
     (XXX explain this in README)
   
   Josh Tobin (tobinrj@tcd.ie), 2015
   ======================================================================== */

// XXX: use Graph class' g6 code instead of repeating it here

#include <iostream>
#include <cmath>
#include <fstream>
#include "graph.hpp"
#include "../clatex/clatex.hpp"
#include "graph_draw.hpp"

using namespace std;
using namespace arma;

unsigned int readN(unsigned long long& r, char* d)
{
  if(d[0] < 126)
  {
    r = (d[0] - 63);
    return 1;
  }
  if(d[1] < 126)
  {
    r = (d[3] - 63) + ((d[2] - 63) << 6) + ((d[1] - 63) << 12);
    return 4;
  }

  r =  (d[4] - 63) + ((d[3] - 63) << 6) + ((d[2] - 63) << 12);
  r = r << 18;
  r += (d[7] - 63) + ((d[6] - 63) << 6) + ((d[5] - 63) << 12);
  return 8;
}

/* Read in the data encoded as *d = R(x), according to
   the graph6 data format specified by McKay.

   b is the preallocated memory to store the resulting data,
   d is the incoming data and n is the length of the incoming
   data in bytes.
*/

void readR(bool* b, char* d, unsigned long long int n)
{
  for(int i=0; i<n; i++)
  {
    char c = d[i];
    c -= 63;
    for(int j=0; j<6; j++)
      b[6*(i+1) - j - 1] = ((c >> j) % 2);
  }
}

double TestProperty(Graph* g)
{
  int n = g->getNumVertices();
  double res = 0.;

  double wk = 0.;
  double neg_sum = 0.;
  double abs_tr_k = 0.;
  double energy = 0.;
  
  for(int i=0; i<n; i++)
  {
    //if(g->spectrum(i) < 0)
    //  res -= pow(g->mainAngle(i),2.) * pow( - g->spectrum(i),3.);
    //else res += 0.5 * pow(g->mainAngle(i),2.) * pow(g->spectrum(i),3.);

    abs_tr_k += pow(fabs(g->spectrum(i)),3.);
    energy += fabs(g->spectrum(i));
    
    if(g->spectrum(i) < 0)
    {
      wk -= pow(g->mainAngle(i),2.) * pow( - g->spectrum(i), 3. );
      neg_sum += pow(g->mainAngle(i),2.) * pow(-g->spectrum(i),3.);
    }
    else
    {
      wk += pow(g->mainAngle(i),2.) * pow(g->spectrum(i), 3. );
    }
  }

  //return 100 - (0.5 * wk - pow(0.5 * energy, 3.));
  //return 100 - (wk - abs_tr_k);
  return energy;
}

double GNP(int n, int trials, double p)
{
  Graph g(n);
  double total = 0;
  for(int j=0; j<trials; j++)
  {
    randomGraphGNP(&g, p, n);
    total += TestProperty(&g);
  }
  total /= trials;
  return total;
}

void GNP_fixp(double* res, int start, int n, int trials, double p)
{
  for(int i=0; i<n; i++)
    res[i] = GNP(start+i, trials, p);
}

void max_graph(Graph* g, int size, string path)
{
  Graph G(1);

  double max_v = 0;
  
  ifstream in;
  in.open(path + to_string(size) + ".g6");
  unsigned long long n;

  char buffer[100];
  bool bytes[1000];
  unsigned int n_read;
  while(!in.eof())
  {
    in.getline(buffer,100);
    if(in.eof())
      break;
    n_read = readN(n, buffer);
 
    G.setNumVertices(n);
    
    unsigned long long nb = n*(n-1) / 2;
    if(nb % 6)
      nb = nb/6 + 1;
    else
      nb = nb/6;
    readR(bytes, buffer+n_read, nb);
 
    unsigned int index = 0;
    for(int i=1; i<n; i++)
    {
      for(int j=0; j<i; j++)
      {
        if(bytes[index++])
          G.addEdge(j,i);
      }
    }

    double tmp = TestProperty(&G);
    if(tmp > max_v)
    {
      max_v = tmp;
      *g = G;
    }
  }
}

void extremal_section(Clatex& report, string title, string file_prefix,
                      int start, int end)
{
  CSection& extSec = report.newSection(title);

  CText* graphCent[end-start+1];
  CDrawing* graphDraw[end-start+1];
  
  for(int i=start; i<=end; i++)
  {
    extSec.addText("$n=" + to_string(i) + "$:\n");
    graphCent[i-start] = & extSec.matchedCmd("center"); 
    graphDraw[i-start] = new CDrawing();
    graphCent[i-start]->addText(graphDraw[i-start]);
  }

  Graph mg(1);
  for(int i=start; i<=end; i++)
  {
    cout << "starting " << i << endl;
    max_graph(&mg, i, file_prefix);
    force_draw(*graphDraw[i-start], mg, 10., 6.);
    cout << "max: " << TestProperty(&mg) << endl;
  }  
}
    
int main()
{

/*-------------------------------------------------
  INITIAL SETUP
  -------------------------------------------------*/
  
  const int num_trials = 500;
  Clatex report;
  report.setTitle("Test Report Document", "grapt");

/*-------------------------------------------------
  SECTION 1:  RANDOM GRAPHS
  -------------------------------------------------*/
  
  CSection& randomSec = report.newSection("Random Graphs");
  randomSec.addText("GNP model with $p=1/2$, $n$ varies:");
  CText& randomCent = randomSec.matchedCmd("center");
  CDrawing* randomDraw = new CDrawing();
  CPlot* randomPlot = new CPlot();
  randomCent.addText(randomDraw);
  randomDraw->addShape(randomPlot);

  int length = 20;
  double data[length];
  GNP_fixp(data, 10, 20, num_trials, 0.5);
  randomPlot->setSize(length);
  for(int i=0; i<20; i++)
  {
    randomPlot->setPoint(i,10+i,data[i]);
  }

  randomSec.addText("GNP model with $n=30$, $p$ varies:");
  CText& randomCent2 = randomSec.matchedCmd("center");
  CDrawing* randomDraw2 = new CDrawing();
  CPlot* randomPlot2 = new CPlot();
  randomCent2.addText(randomDraw2);
  randomDraw2->addShape(randomPlot2);

  randomPlot2->setSize(1./0.05);
  for(int i=0; i < 1. / 0.05; i++)
  {
    randomPlot2->setPoint(i,i*0.05,GNP(30,num_trials,i*0.05));
  }

/*-------------------------------------------------
  SECTION 2:  EXTREMAL GRAPHS
  -------------------------------------------------*/
  
  extremal_section(report, "Extremal Connected Graphs", "connected/con", 4, 8);    
  extremal_section(report, "Extremal Trees", "trees/tree", 4, 18);  
  extremal_section(report, "Regular Graphs", "reg/reg", 2, 12);   

  
/*-------------------------------------------------
  GENERATE REPORT
  -------------------------------------------------*/
  
  report.write("rep.tex");
  return 0;
}

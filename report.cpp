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

void ext_graph(Graph& min_g, Graph& max_g, double& min_val,
               double& max_val, double& avg_val, int size, string path)
{
  Graph G(1);

  bool first_graph = true;
  
  ifstream in;
  in.open(path + to_string(size) + ".g6");
  unsigned long long n;

  char buffer[100];
  bool bytes[1000];
  unsigned int n_read;
  unsigned int num_graphs = 0;
  avg_val = 0;
  
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
    avg_val += tmp;
    num_graphs++;
    if(first_graph || tmp < min_val)
    {
      min_val = tmp;
      min_g = G;
    }
    if(first_graph || tmp > max_val)
    {
      max_val = tmp;
      max_g = G;
    }
    if(first_graph)
    {
      first_graph = false;
    }
  }

  avg_val /= num_graphs;
}

void extremal_section(Clatex& report, string title, string file_prefix,
                      int start, int end)
{
  CSection& extSec = report.newSection(title);

  CDrawing* minDraw[end-start+1];  
  CDrawing* maxDraw[end-start+1];
  CText* avgText[end-start+1];
  CText* minText[end-start+1];
  CText* maxText[end-start+1];

  CText* avgVal[end-start+1];
  CText* minVal[end-start+1];
  CText* maxVal[end-start+1];

  
  for(int i=start; i<=end; i++)
  {
    vector<CText*> row[3];
    row[0].resize(3);
    row[1].resize(3);
    row[2].resize(3);

    minDraw[i-start] = new CDrawing();
    maxDraw[i-start] = new CDrawing();

    
    minText[i-start] = new CText();
    maxText[i-start] = new CText();
    avgText[i-start] = new CText();

    minVal[i-start] = new CText();
    maxVal[i-start] = new CText();
    avgVal[i-start] = new CText();    
    
    extSec.addText("$n=" + to_string(i) + "$:\\\\ \\\\ \n");
    //extSec.addText(avgText[i-start]);

    CTable* table = new CTable();
    table->setNumCols(3);
    table->setJust("m{2cm} m{4cm} m{5cm}");
    extSec.addText("\\fbox{"); // XXX abstract this to Clatex
    extSec.addText(table);
    extSec.addText("}\\\\ \\\\");

    row[0][0] = avgText[i-start];
    row[0][1] = avgVal[i-start];
    row[0][2] = new CText("");
    table->addRow(row[0]);
        
    row[1][0] = minText[i-start];
    row[1][1] = minVal[i-start];
    row[1][2] = minDraw[i-start];
    table->addRow(row[1]);

    row[2][0] = maxText[i-start];
    row[2][1] = maxVal[i-start];
    row[2][2] = maxDraw[i-start];
    table->addRow(row[2]);
    
//    extSec.addText(minText[i-start]);
//    CText& minCent = extSec.matchedCmd("center"); 
//    minCent.addText(minDraw[i-start]);

//    extSec.addText(maxText[i-start]);
//    CText& maxCent = extSec.matchedCmd("center"); 
//    maxDraw[i-start] = new CDrawing();
//    maxCent.addText(maxDraw[i-start]);
  }

  Graph maxG(1), minG(1);
  double max_val, min_val, avg_val;
  for(int i=start; i<=end; i++)
  {
    cout << "starting " << i << endl;
    ext_graph(minG, maxG, min_val, max_val, avg_val, i, file_prefix);
    force_draw(*minDraw[i-start], minG, 3., 3., 2000);
    force_draw(*maxDraw[i-start], maxG, 3., 3., 2000);
    avgText[i-start]->addText("Avg value:");
    avgVal[i-start]->addText(to_string(avg_val));
    minText[i-start]->addText("Min value:");
    minVal[i-start]->addText(to_string(min_val));
    maxText[i-start]->addText("Max value:");
    maxVal[i-start]->addText(to_string(max_val));
    cout << "min: " << min_val << " max: " << max_val << " avg: " << avg_val << endl;
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
  
  extremal_section(report, "Extremal Connected Graphs", "connected/con", 4, 7);    
  extremal_section(report, "Extremal Trees", "trees/tree", 4, 15);  
  extremal_section(report, "Regular Graphs", "reg/reg", 2, 10);   

  
/*-------------------------------------------------
  GENERATE REPORT
  -------------------------------------------------*/
  
  report.write("rep.tex");
  return 0;
}

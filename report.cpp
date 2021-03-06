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

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>

#include <pthread.h>

#include "graph.hpp"
#include "../clatex/clatex.hpp"
#include "graph_draw.hpp"

using namespace std;
using namespace arma;

double TestProperty(Graph* g)
{
  int n = g->getNumVertices();
  return 1. - pow(g->mainAngle(n-1),2.);
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

struct ext_graph_args {
  ext_graph_args(Graph& _min_g, Graph& _max_g, double& _min_v,
                 double& _max_v, double& _avg_v, string& _graphs)
   : min_g(_min_g), max_g(_max_g), min_val(_min_v), max_val(_max_v), avg_val(_avg_v), graphs(_graphs) {};
  Graph& min_g;
  Graph& max_g;
  double& min_val;
  double& max_val;
  double& avg_val;
  string& graphs;
};

void* ext_graph_slave(void* _args)
{
  ext_graph_args* args = (ext_graph_args*) _args;
  Graph& min_g = args->min_g;
  Graph& max_g = args->max_g;
  double& min_val = args->min_val;
  double& max_val = args->max_val;
  double& avg_val = args->avg_val;
  
  bool first_graph = true;

  
  istringstream in(args->graphs);

  unsigned long long n;
  Graph G;
  
  char buffer[100];
  bool bytes[1000];
  unsigned int n_read;
  unsigned int num_graphs = 0;
  avg_val = 0;
  
  while(in)
  {
    in.getline(buffer,100);
    if(strlen(buffer) == 0 && !in)
      break;
    G.from_g6(buffer);
 
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
  if(num_graphs > 0)
    avg_val /= num_graphs;
  return 0;
}

void ext_graph(Graph& min_g, Graph& max_g, double& min_val,
               double& max_val, double& avg_val, int size, string path)
{
  const int n_threads = 3;
  string s[n_threads];

  string line;
  
  ifstream in;
  string file_path = path + to_string(size) + ".g6";
  in.open(file_path); 
  if(in.fail())
  {
    cerr << "Failed to open " << file_path << endl;
    return;
  }  
  int index = 0;
  while(!in.eof())
  {
    in >> line;
    if(line.length() == 0 || !in)
      break;
    s[index % n_threads] += line; // CURP
    s[index % n_threads] += "\n";
    index++;
  }

  Graph tmin_g[n_threads];
  Graph tmax_g[n_threads];  
  double tmin_v[n_threads];
  double tmax_v[n_threads];
  double tavg_v[n_threads]; 


  
  pthread_t thread[n_threads];
  ext_graph_args* args[n_threads];
  for(int i=0; i<n_threads; i++)
  {
    args[i] = new ext_graph_args(tmin_g[i], tmax_g[i], tmin_v[i], tmax_v[i], tavg_v[i], s[i]);
    pthread_create(&thread[i], NULL, ext_graph_slave, (void*) args[i]);
    //ext_graph_slave((void*) args[i]);
  }

  for(int i=0; i<n_threads; i++)
  {
    pthread_join(thread[i], NULL);
  }
  min_val = tmin_v[0];
  max_val = tmax_v[0];
  avg_val = 0;
  
  for(int i=0; i<n_threads; i++)
  {
    if(tmin_v[i] < min_val)
    {
      min_val = tmin_v[i];
      min_g = tmin_g[i];
    }
    if(tmax_v[i] > max_val)
    {
      max_val = tmax_v[i];
      max_g = tmax_g[i];
    }

    avg_val += tavg_v[i] / n_threads;
    
    delete args[i];
  }
}

void extremal_section(Clatex& report, string title, string file_prefix,
                      int start, int end)
{
  CSection& extSec = report.newSection(title, true);

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
  }

  Graph maxG(1), minG(1);
  double max_val, min_val, avg_val;

  cout << "Starting data set " << file_prefix << endl;
  
  for(int i=start; i<=end; i++)
  {
    cout << '\t' << i;
    ext_graph(minG, maxG, min_val, max_val, avg_val, i, file_prefix);
    force_draw(*minDraw[i-start], minG, 3., 3., 4000);
    force_draw(*maxDraw[i-start], maxG, 3., 3., 4000);
    avgText[i-start]->addText("Avg value:");
    avgVal[i-start]->addText(to_string(avg_val));
    minText[i-start]->addText("Min value:");
    minVal[i-start]->addText(to_string(min_val));
    maxText[i-start]->addText("Max value:");
    maxVal[i-start]->addText(to_string(max_val));
    cout << ": min=" << min_val << ", max=" << max_val << ", avg=" << avg_val << endl;
  }  
}
    
int main()
{

/*-------------------------------------------------
  INITIAL SETUP
  -------------------------------------------------*/
  
  const int num_trials = 20; // 2000
  Clatex report;
  report.setTitle("Test Report Document", "grapt");
  report.generateTOC(true);
  
/*-------------------------------------------------
  SECTION 1:  RANDOM GRAPHS
  -------------------------------------------------*/
  
  CSection& randomSec = report.newSection("Random Graphs", true);
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
  
  extremal_section(report, "Extremal Connected Graphs", "connected/con", 4, 9); // max: 10    
  extremal_section(report, "Extremal Trees", "trees/tree", 4, 19);              // max: 22  
  extremal_section(report, "Regular Graphs", "reg/reg", 2, 12);                  // max: 14 
  extremal_section(report, "$C_4$--free Graphs", "c4free/cf", 4, 12);           // max: 15   
  extremal_section(report, "Bipartite Graphs", "bipartite/bip", 4, 12);         // max: 14

  
/*-------------------------------------------------
  GENERATE REPORT
  -------------------------------------------------*/
  
  report.write("rep.tex");
  return 0;
}

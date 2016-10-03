/* ========================================================================
   grapt:
     Generates a LaTeX report of a user-defined graph property.

   graph_draw.cpp:  Force-directed graph drawing functions. 
   
   Josh Tobin (tobinrj@tcd.ie), 2015
   ======================================================================== */

#include "graph_draw.hpp"

#define EPS 1e-9

static void out128(__m128 t)
{
  float* f = (float*) &t;
  cout << *f << ' ';
  f++;
  cout << *f << ' ';
  f++;
  cout << *f << ' ';
  f++;
  cout << *f << endl;
}

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
void sse_force_draw(CDrawing& drawing, Graph& g, double x_size, double y_size, int num_iterations)
{
  int n = g.getNumVertices();
  int N = n + ((4 - n%4) % 4);

  float c1 = 2., c2 = 1., c3 = 1., c4 = .05;

  float _neg_c1 = -c1;
  float _inv_c2 = 1. / c2;
  float small_cons = 0.01;

  __m128 neg_c1 = _mm_load_ps1(&_neg_c1);  
  __m128 inv_c2 = _mm_load_ps1(&_inv_c2);
  __m128 m_c3 = _mm_load_ps1(&c3);
  __m128 m_c4 = _mm_load_ps1(&c4);
  __m128 m_sc = _mm_load_ps1(&small_cons);
  
  float * __attribute__ ((aligned (16))) pts_x = new float[N] ;  
  float * __attribute__ ((aligned (16))) pts_y = new float[N];
  float * __attribute__ ((aligned (16))) forces_x = new float[N] ;  
  float * __attribute__ ((aligned (16))) forces_y = new float[N];

//  float* dist  = new float[N] __attribute__ ((aligned (16)));

  float __attribute__ ((aligned (16))) *adj[n];
  float __attribute__ ((aligned (16))) *inv_adj[n];
  
  for(int i=0; i<n; i++)
  {
    adj[i] = new float[N];
    inv_adj[i] = new float[N];
    
    for(int j=0; j<n; j++)
    {
      adj[i][j] = (g.isConnected(i,j) ? 1. : 0.);
      inv_adj[i][j] = 1. - adj[i][j];      
    }
    for(int j=n; j<N; j++)
    {
      adj[i][j] = 0.;
      inv_adj[i][j] = 0.;
    }
  }

  __m128 tmp1;
  __m128 tmp2;
  __m128 tmp3;
  __m128 tmp4;
 
  for(int i=0; i<n; i++)
  {
    pts_x[i] = drand48() * (x_size / 1.2);
    pts_y[i] = drand48() * (y_size / 1.2);
  }

  for(int i=0; i<num_iterations; i++)
  {
    for(int j=0; j<n; j++)
    {
      __m128* ppts_x   = (__m128*) pts_x;
      __m128* ppts_y   = (__m128*) pts_y;
      __m128* padj     = (__m128*) adj[j];
      __m128* pinv_adj = (__m128*) inv_adj[j];


      tmp1 = _mm_load_ps1(pts_x + j);
      tmp2 = _mm_load_ps1(pts_y + j);

      __m128 fc, fc1, fc2, dist, diff_x, diff_y;

      float sum_x = 0., sum_y = 0.;    

      forces_x[j] = forces_y[j] = 0;
      
      while(ppts_x < ((__m128*) pts_x) + N/4)
      {
        tmp3 = _mm_sub_ps(*ppts_x, tmp1);
        dist = _mm_mul_ps(tmp3, tmp3);

        tmp3 = _mm_sub_ps(*ppts_y, tmp2);
        dist = _mm_add_ps(dist, _mm_mul_ps(tmp3, tmp3));

        dist = _mm_sqrt_ps(dist);
        dist = _mm_add_ps(dist, m_sc);

        fc1 = _mm_mul_ps(dist, inv_c2);

        
        float* tmp_f = (float*) &fc1;
        *tmp_f = log(*tmp_f); tmp_f++;
        *tmp_f = log(*tmp_f); tmp_f++;
        *tmp_f = log(*tmp_f); tmp_f++;
        *tmp_f = log(*tmp_f); 
 
        fc1 = _mm_mul_ps(*padj, fc1);
        
        fc1 = _mm_mul_ps(fc1, neg_c1);

        fc2 = _mm_rcp_ps(dist);

        fc2 = _mm_mul_ps(*pinv_adj, fc2);
        fc2 = _mm_mul_ps(fc2, fc2);
        fc2 = _mm_mul_ps(fc2, m_c3);

 
        
        fc = _mm_add_ps(fc1, fc2);
        fc = _mm_mul_ps(fc, m_c4);

        diff_x = _mm_sub_ps(tmp1, *ppts_x);
        diff_x = _mm_div_ps(diff_x, dist);
        diff_x = _mm_mul_ps(diff_x, fc);
        diff_y = _mm_sub_ps(tmp2, *ppts_y);
        diff_y = _mm_div_ps(diff_y, dist);
        diff_y = _mm_mul_ps(diff_y, fc);


        tmp_f = (float*) &diff_x;
        sum_x += *tmp_f; tmp_f++;
        sum_x += *tmp_f; tmp_f++;
        sum_x += *tmp_f; tmp_f++;
        sum_x += *tmp_f;
        
        tmp_f = (float*) &diff_y;
        sum_y += *tmp_f; tmp_f++;
        sum_y += *tmp_f; tmp_f++;
        sum_y += *tmp_f; tmp_f++;
        sum_y += *tmp_f; 

        sum_x += *((float*) &diff_x);
        sum_y += *((float*) &diff_y);
                
        ppts_x++;
        ppts_y++;
        padj++;
        pinv_adj++;
      }

      forces_x[j] += sum_x;
      forces_y[j] += sum_y;
    }

    __m128* ppts_x   = (__m128*) pts_x;
    __m128* ppts_y   = (__m128*) pts_y;
    __m128* pfcs_x   = (__m128*) forces_x;
    __m128* pfcs_y   = (__m128*) forces_y;

    while(ppts_x < ((__m128*) pts_x) + N/4)
    {
      _mm_store_ps((float*)ppts_x, _mm_add_ps(*ppts_x, *pfcs_x));
      _mm_store_ps((float*)ppts_y, _mm_add_ps(*ppts_y, *pfcs_y));
      
      ppts_x++;
      ppts_y++;
      pfcs_x++;
      pfcs_y++;
    }    
  }

  CPoint* pts[n];
  for(int i=0; i<n; i++)
  {
    pts[i] = new CPoint(pts_x[i], pts_y[i]);
  }
  
  double max_x = pts_x[0], min_x = pts_x[0];
  double max_y = pts_y[0], min_y = pts_y[0];
  
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
  if(width < x_size && height < y_size)
    width = x_size = height = y_size = 1.;

  double scale_x = x_size / width;
  double scale_y = y_size / height;
  if(scale_x > scale_y)
    scale_x = scale_y;
  else
    scale_y = scale_x;
  
  
  for(int i=0; i<n; i++)
  {
    pts[i]->x = (pts[i]->x - min_x) * scale_x;
    pts[i]->y = (pts[i]->y - min_y) * scale_y;    
  }
  
  for(int i=0; i<n; i++)
  {
    drawing.addShape(pts[i]);
  }
  
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
      if(g.isConnected(i,j))
        drawing.addShape(new CLine(*pts[i],*pts[j]));  

  delete[] pts_x;
  delete[] pts_y;
  delete[] forces_x;
  delete[] forces_y;

  for(int i=0; i<n; i++)
  {
    delete[] adj[i];
    delete[] inv_adj[i];
  }
    
}

void force_draw(CDrawing& drawing, Graph& g, double x_size, double y_size, int num_iterations, std::string* labels)
{
  int n = g.getNumVertices();

  float c1 = 2., c2 = 1., c3 = 1., c4 = .05;  
  CPoint* pts[n];
  CPoint* forces[n];
  for(int i=0; i<n; i++)
  {
    forces[i] = new CPoint(0,0);
    pts[i] = new CPoint(drand48() * (double) (x_size / 1.2), drand48() * (double) (y_size / 1.2));
  }

  for(int i=0; i<num_iterations; i++)
  {
    for(int j=0; j<n; j++)
    {
      double fx = 0.,  fy = 0., dist, fc;
      forces[j]->x = forces[j]->y = 0;
      for(int k=0; k<n; k++)
      {
        //if(j==k)
        //  continue;
        dist = 0.01 + sqrt(pow(pts[j]->x - pts[k]->x,2.) + pow(pts[j]->y - pts[k]->y,2.));
        if(g.isConnected(j,k))
          fc = - c1 * log(dist / c2);
        else
          fc = c3 / pow(dist,2.);
        
        fc *= c4;        
        
        fx += fc * (pts[j]->x - pts[k]->x) / dist;
        fy += fc * (pts[j]->y - pts[k]->y) / dist;
      }

      forces[j]->x += fx;
      forces[j]->y += fy;
    }
    for(int j=0; j<n; j++)
    {
      pts[j]->x += forces[j]->x;
      pts[j]->y += forces[j]->y;      
    }
  }

  double max_x = pts[0]->x, min_x = pts[0]->x;
  double max_y = pts[0]->x, min_y = pts[0]->x;
  
  for(int i=0; i<n; i++)
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
  if(width < x_size && height < y_size)
    width = x_size = height = y_size = 1.;

  double scale_x = x_size / width;
  double scale_y = y_size / height;
  if(scale_x > scale_y)
    scale_x = scale_y;
  else
    scale_y = scale_x;
  
  
  for(int i=0; i<n; i++)
  {
    pts[i]->x = (pts[i]->x - min_x) * scale_x;
    pts[i]->y = (pts[i]->y - min_y) * scale_y;    
  }
  
  for(int i=0; i<n; i++)
  {
    if(!labels)
      drawing.addShape(pts[i]);
    else
    {
      CLabel* label = new CLabel();
      label->x = pts[i]->x;
      label->y = pts[i]->y;
      label->text = labels[i];
      label->name = "v" + to_string(i);
      //delete pts[i];
      drawing.addShape(label);
    }
  }
  
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
      if(g.isConnected(i,j))
      {
        if(!labels)
          drawing.addShape(new CLine(*pts[i],*pts[j]));
        if(labels)
          drawing.addShape(new CLine("v" + to_string(i), "v" + to_string(j)));
      }

  if(labels)
    for(int i=0; i<n; i++)
      delete pts[i];
  
  for(int i=0; i<n; i++)
    delete forces[i];
}

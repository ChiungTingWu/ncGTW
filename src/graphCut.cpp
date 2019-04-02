#include <Rcpp.h>
#include "ibfs.h"
using namespace Rcpp;


double round(double a)
{
  return floor(a + 0.5);
}


int isInteger(double a)
{
  return (abs(a - round(a)) < 1e-6);
}

typedef IBFSGraph<double, double, double> GraphType;

//[[Rcpp::export]]
NumericVector graphCut(NumericMatrix termW, NumericMatrix edges)
{
  //node and edge number
  int numNodes = termW.nrow();
  int numEdges = edges.nrow();

  //prepare graph
  GraphType *g = new GraphType( numNodes, numEdges);
  for(int i = 0; i < numNodes; ++i)
  {
    g -> add_node(1);
    g -> add_tweights( i, termW[i], termW[numNodes + i]);
  }

  for(int i = 0; i < numEdges; ++i)
  {
    if(edges[i] < 1 || edges[i] > numNodes || edges[numEdges + i] < 1 || edges[numEdges + i] > numNodes || edges[i] == edges[numEdges + i] || !isInteger(edges[i]) || !isInteger(edges[numEdges + i])){
      stop("graphCut:pairwisePotentials", "Some edge has invalid vertex numbers and therefore it is ignored");
    }
    else{
      if(edges[2 * numEdges + i] + edges[3 * numEdges + i] < 0){
        stop("graphCutMex:pairwisePotentials", "Some edge is non-submodular and therefore it is ignored");
      }
      else
      {
        if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] >= 0)
          g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), edges[2 * numEdges + i], edges[3 * numEdges + i]);
        else
          if (edges[2 * numEdges + i] <= 0 && edges[3 * numEdges + i] >= 0)
          {
            g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i] + edges[2 * numEdges + i]);
            g -> add_tweights((GraphType::node_id)round(edges[i] - 1), 0, edges[2 * numEdges + i]);
            g -> add_tweights((GraphType::node_id)round(edges[numEdges + i] - 1),0 , -edges[2 * numEdges + i]);
          }
          else
            if (edges[2 * numEdges + i] >= 0 && edges[3 * numEdges + i] <= 0)
            {
              g -> add_edge((GraphType::node_id)round(edges[i] - 1), (GraphType::node_id)round(edges[numEdges + i] - 1), edges[3 * numEdges + i] + edges[2 * numEdges + i], 0);
              g -> add_tweights((GraphType::node_id)round(edges[i] - 1),0 , -edges[3 * numEdges + i]);
              g -> add_tweights((GraphType::node_id)round(edges[numEdges + i] - 1), 0, edges[3 * numEdges + i]);
            }
            else{
              stop("graphCut:pairwisePotentials", "Something strange with an edge and therefore it is ignored");
            }
      }
    }
  }
  //compute flow
  double flow = g -> maxflow();

  NumericVector segment(numNodes + 1);

  //output minimum cut
  for(int i = 0; i < numNodes; ++i)
    segment[i] = g -> what_segment(i);

  //output minimum value
  segment[numNodes] = flow;

  delete g;

  return segment;
}



#include <RcppArmadillo.h>
#include <iostream>
#include <tictoc.h>

using namespace arma;
 

double intermediation_cpp(const mat& adj,const mat& radj) {
  mat nradj=adj-radj; //nonreciprocal
  vec giving = sum(nradj,1);
  //std::cout << giving <<endl;
  vec receiving = trans(sum(nradj,0));
  //std::cout << receiving <<endl;
  vec giving_and_receiving = min(giving,receiving);
  vec giving_or_receiving = max(giving,receiving);
  vec giving_and_receiving_ =  nonzeros(giving_and_receiving);
  vec giving_or_receiving_ =  nonzeros(giving_or_receiving);
  //std::cout <<"a"<< giving_and_receiving_ <<endl;
  //std::cout <<"av"<< giving_and_receiving_.n_elem <<endl;
  //std::cout <<"aa"<< ((float)giving_and_receiving_.n_elem)/((float)giving_or_receiving_.n_elem) << endl;
  if (giving_or_receiving_.n_elem==0) return(0);
  return(((float)giving_and_receiving_.n_elem)/((float)giving_or_receiving_.n_elem));
}
// [[Rcpp::export]]
double recip_cpp(const mat& adj, const mat& radj) {
  return(accu(radj)/accu(adj));
}

double support_fast2_cpp(const mat& m,const mat& undir) {
  uvec linksifsupported = find(undir>0);
  if (linksifsupported.n_elem==0) return(0);
  mat twopaths = undir*undir;
  return(mean(vectorise(clamp(twopaths.elem(linksifsupported),0,1))));
}



void breadthfirst(uword i,uword & c, uvec & membership, const mat & adj) {
  uvec sc=find(adj.row(i)); //the search space could be rstricted further
  for(uword j : sc) {
    if (membership(j)!=0) {continue;}
    membership(j)=c;
    breadthfirst(j, c, membership, adj);
  }
  return;
}
// [[Rcpp::export]]
double component_counts(const mat& transferstructure) {  
  mat adj = transferstructure+trans(transferstructure);
  uword t_components_no=0;
  uvec t_components_membership = zeros<uvec>(transferstructure.n_rows);
  for(uword i=0; i<transferstructure.n_rows; i++) {
    if (t_components_membership(i)!=0) {continue;} else {t_components_membership(i)=++t_components_no;}
    breadthfirst(i,t_components_no,t_components_membership,adj);
  }
  return t_components_no; 
}

// [[Rcpp::export]]
double Ccomponents(const mat& transferstructure,uvec & t_components_csize,uvec & t_components_membership, mat & t_components_matrix, mat & t_conmat) {  
  mat adj = transferstructure+trans(transferstructure);
  uword t_components_no=0;
  t_components_membership = zeros<uvec>(transferstructure.n_rows);
  for(uword i=0; i<transferstructure.n_rows; i++) {
    if (t_components_membership(i)!=0) {continue;} else {t_components_membership(i)=++t_components_no;}
    breadthfirst(i,t_components_no,t_components_membership,adj);
  }
  t_components_matrix=zeros(transferstructure.n_cols,t_components_no);
  t_conmat=zeros(transferstructure.n_rows,transferstructure.n_cols);
  t_components_csize=uvec(t_components_no);
  for(uword cc=1; cc<=t_components_no;cc++) {
    t_conmat(find(t_components_membership==cc),find(t_components_membership==cc)).ones();
    uvec zz = find(t_components_membership==cc);
    uvec aa = {cc-1};
    t_components_matrix.submat(zz,aa).ones();
    t_components_csize(cc-1)=zz.n_elem;
  }    
  
  return t_components_no;
} 
// [[Rcpp::export]]

double forestness_cpp(const mat& adj) {
  return(2*(adj.n_rows-component_counts(adj)) / accu(adj));
}



//shorest_path_code from here: https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/

// The program is for adjacency matrix representation of the graph 

// A utility function to find the vertex with minimum distance value, from 
// the set of vertices not yet included in shortest path tree 
int minDistance(vec dist, vec sptSet) { 
  // Initialize min value 
  double min = datum::inf, min_index = 0; 
  int V=dist.n_elem;

  for (int v = 0; v < V; v++) {

    if (sptSet(v) == 0 && dist(v)<=min) 
      min = dist(v), min_index = v; 
  }
  return min_index; 
} 



//Breath first shortest distances: https://www.geeksforgeeks.org/multi-source-shortest-path-in-unweighted-graph/
// This function returns for a given adjacency matrix and a given source node the shortest distances to all other nodes
// [[Rcpp::export]]
mat zeroOneBFS(const mat& m,int src) 
{ 
  src-=1;
  uword V = m.n_rows;
  // Initialize distances from given source 
  vec dist(V);
  dist.fill(datum::inf);

  // double ende queue to do BFS. 
  std::deque <int> Q; 
  dist(src) = 0; 
  Q.push_back(src); 
  
  while (!Q.empty()) 
  { 
    int v = Q.front(); 
    Q.pop_front(); 
    rowvec vsrow=m.row(v);
    uvec edges = find(vsrow);
    for (uword i=0; i<edges.n_elem; i++) 
    { 
      // checking for the optimal distance 
      if (dist(edges(i)) > dist(v) + m(v,edges(i))) 
      { 
        dist(edges(i)) = dist(v) + m(v,edges(i)); 
        
        // Put 0 weight edges to front and 1 weight 
        // edges to back so that vertices are processed 
        // in increasing order of weights. 
        if (m(v,edges(i)) == 0) 
          Q.push_front(edges(i)); 
        else
          Q.push_back(edges(i)); 
      } 
    } 
  } 
  
  return(dist);
} 


// [[Rcpp::export]]
mat BFS_dist_all(mat graph) 
{ 
  int V=graph.n_rows;
  mat ret = zeros(V,V);
  for (int i=1;i<=V;i++)
    ret.row(i-1)=trans(zeroOneBFS(graph,i));
  return(ret);
}







// Function that implements Dijkstra's single source shortest path algorithm 
// for a graph represented using adjacency matrix representation 
// [[Rcpp::export]]
vec dijkstra(mat graph, int src) 
{ 
  src-=1;
  int V=graph.n_rows;
  graph.replace(0,V+1);
  vec dist=zeros(V); // The output array.  dist[i] will hold the shortest 
  // distance from src to i 
  
  vec sptSet=zeros(V); // sptSet[i] will be true if vertex i is included in shortest 
  // path tree or shortest distance from src to i is finalized 
  
  // Initialize all distances as INFINITE and stpSet[] as false 
  for (int i = 0; i < V; i++) 
    dist(i) = datum::inf, sptSet(i) = false; 
  
  // Distance of source vertex from itself is always 0 
  dist(src) = 0; 
  
  // Find shortest path for all vertices 
  for (int count = 0; count < V - 1; count++) { 
    // Pick the minimum distance vertex from the set of vertices not 
    // yet processed. u is always equal to src in the first iteration. 
    int u = minDistance(dist, sptSet); 

    // Mark the picked vertex as processed 
    sptSet(u) = 1;

    // Update dist value of the adjacent vertices of the picked vertex. 
    for (int v = 0; v < V; v++)  {

      // Update dist[v] only if is not in sptSet, there is an edge from 
      // u to v, and total weight of path from src to  v through u is 
      // smaller than current value of dist[v] 
      if (sptSet(v)==0 && graph(u,v)!=0 && dist(u) != datum::inf 
            && dist(u) + graph(u,v) < dist(v)) 
        dist(v) = dist(u) + graph(u,v); 
    }
  } 
  dist.replace(V+1,datum::inf);
  
  // print the constructed distance array 
  return(dist); 
} 

// [[Rcpp::export]]
mat dijkstra_all(mat graph) 
{ 
  int V=graph.n_rows;
  mat ret = zeros(V,V);
  for (int i=1;i<=V;i++)
    ret.row(i-1)=trans(dijkstra(graph,i));
  return(ret);
}

// [[Rcpp::export]]
double dskewness(const mat& adj)
{
  
  vec degree_distribution = sum(adj,1);
  return(mean(pow((degree_distribution-mean(degree_distribution))/stddev(degree_distribution),3)));
  
}
// [[Rcpp::export]]
double degreesd(const mat& adj)
{
  
  vec degree_distribution = sum(adj,1);
  return(stddev(degree_distribution));
  
}

// [[Rcpp::export]] 
vec compute_moments_cpp(const mat& btransfers,const mat& kinship,const mat& distance,const vec& income) {
  vec values=btransfers.elem(find(btransfers>0));
  if (values.n_elem>0) {
      if ((values.min()<1) | (values.max()>1)) {
      Rcpp::Rcout << "your matrix is not binary";
    }
  }
  uvec offdiag=find(eye(kinship.n_rows,kinship.n_rows)==0);
  
  mat m_undir=max(btransfers,trans(btransfers));
  mat m_recip=btransfers%trans(btransfers);
  
  double fb2=forestness_cpp(m_undir);
  double ib=intermediation_cpp(btransfers,m_recip);
  
  
  mat con = (1/income)*trans(income);
  mat logcon = log(trans(con));

  /*resulst in the current draft are based on this, but for the current run we don't need to compute these
   * 
  mat dat2(offdiag.n_elem,2);
  dat2.col(0)=abs(vectorise(logcon(offdiag)));
  dat2.col(1)=vectorise(kinship(offdiag));
  int n2 = dat2.n_rows, k2 = dat2.n_cols;
  //lm(pmax(transfers,t(transfers))[offdiag] ~ -1 + kinship[offdiag] + abs(logcon[offdiag]))
  vec coef2 = solve(dat2, m_undir(offdiag)); 
  vec resid2 = m_undir(offdiag) - dat2*coef2; 
  double sig22 = as_scalar(trans(resid2)*resid2/(n2-k2));
  */
  //auxiliary regression for all
  mat dat_full(offdiag.n_elem,3);
  dat_full.col(0)=ones(offdiag.n_elem);
  dat_full.col(1)=vectorise(logcon(offdiag));
  dat_full.col(2)=vectorise(kinship(offdiag));
  
  //int n_full = dat_full.n_rows, k_full = dat_full.n_cols;
  vec coef_full = solve(dat_full, btransfers(offdiag)); 
  //vec resid_full = btransfers(offdiag) - dat_full*coef_full; 
  //double sig_full = as_scalar(trans(resid_full)*resid_full/(n_full-k_full));
  
  //double degree_skewness = dskewness(m_undir);
  double degree_sd = degreesd(m_undir);
  /*
   vec degree_distribution = sum(m_undir,1);
  */
  
  //create datinteract directly from dat_Full with an additioncal column with the interaction of col 2 and 3
  mat datinteract(dat_full.n_rows,dat_full.n_cols+1);
  datinteract.cols(0,dat_full.n_cols-1)=dat_full;
  datinteract.col(dat_full.n_cols)=dat_full.col(1)%dat_full.col(2);
  vec coefinteract = solve(datinteract, btransfers(offdiag));
  
  /* correlation between having multiple paths and distance
  mat twopaths=(m_undir*m_undir)%m_undir;
  mat corrdoubledist_=cor(twopaths(offdiag),distance(offdiag));
  double corrdoubledista=corrdoubledist_(0);
   */
  
  /*correlation between distance-degree (centrality) and link degree: The idea: s/o with few relatives, will have high constraints and thus more links?
  vec dist_distribution = sum(distance,1);
  vec outdegree_distribution = sum(btransfers,1);
  mat correlation_of_degrees = cor(degree_distribution,dist_distribution);
  */
//  1  2  3  7  10 13
  vec ret = {degree_sd, //density, //1 
             fb2, //2
             ib,  //3
             coef_full(0),//4
             coef_full(1),//5,
             coef_full(2),//6
             99,//c1(0),//7
             99, //degree_skewness, //c2(0),//8
             99, //coef2(0),//9
             coefinteract(0),//sig22,//10
             coefinteract(1),//coef2(1),//correlation_of_degrees(0), //// 11
             coefinteract(2), //12
             coefinteract(3)//degree_skewness //13
    };
  return(ret.replace(datum::nan,0));
}
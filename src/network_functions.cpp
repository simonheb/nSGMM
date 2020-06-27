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

double forestness_cpp(const mat& adj) {
  return(2*(adj.n_rows-component_counts(adj)) / accu(adj));
}



//shorest_path_code from here: https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/

// The program is for adjacency matrix representation of the graph 

// A utility function to find the vertex with minimum distance value, from 
// the set of vertices not yet included in shortest path tree 
int minDistance(vec dist, vec sptSet) { 
  // Initialize min value 
  double min = datum::inf, min_index; 
  int V=dist.n_elem;

  for (int v = 0; v < V; v++) {

    if (sptSet(v) == 0 && dist(v)<=min) 
      min = dist(v), min_index = v; 
  }
  return min_index; 
} 



//Breath first shortest distances: https://www.geeksforgeeks.org/multi-source-shortest-path-in-unweighted-graph/
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
vec compute_moments_cpp(const mat& btransfers,const mat& kinship,const mat& distance,const vec& income,const vec& theta) {
  vec values=btransfers.elem(find(btransfers>0));
  if (values.n_elem>0) {
    if (values.min()<1 | values.max()>1) {
      Rcpp::Rcout << "your matrix is not binary";
    }
  }
  uvec offdiag=find(eye(kinship.n_rows,kinship.n_rows)==0);
  
  mat m_undir=max(btransfers,trans(btransfers));
  mat m_recip=btransfers%trans(btransfers);
  mat pl=BFS_dist_all(m_undir);
  pl=pl.replace(0,datum::nan);
  pl=pl.replace(datum::inf,btransfers.n_rows);
  double pathlenghts=mean(pl.elem(find_finite(pl)))/btransfers.n_rows;
  double fb2=forestness_cpp(m_undir);
  double ib=intermediation_cpp(btransfers,m_recip);
  double sa=support_fast2_cpp(btransfers,m_undir);
  double ra=recip_cpp(m_undir,m_recip);
  mat c1=cor(vectorise(m_undir(offdiag)),vectorise(kinship(offdiag)));
  mat c2=cor(vectorise(m_undir(offdiag)),vectorise(distance(offdiag)));
  double density=mean(vectorise(btransfers(offdiag)));
  
  mat con = (1/income)*trans(income);
  mat logcon = log(trans(con));
  double c3;
  if (density==0) {
    c3=0;
  }else{
    c3=mean(abs(logcon(find(m_undir))));
  }
  mat dat(offdiag.n_elem,7);
  dat.col(0)=vectorise(logcon(offdiag));
  dat.col(1)=sign(dat.col(0));
  dat.col(2)=vectorise(kinship(offdiag));
  dat.col(3)=dat.col(0)%dat.col(0);
  dat.col(4)=dat.col(0)%dat.col(1);
  dat.col(5)=dat.col(0)%dat.col(2);
  dat.col(6)=dat.col(1)%dat.col(2);

  int n = dat.n_rows, k = dat.n_cols;
  
  vec coef = solve(dat, btransfers(offdiag)); 
  vec resid = btransfers(offdiag) - dat*coef; 
  
  double sig2 = as_scalar(trans(resid)*resid/(n-k));
  
  
  mat dat2(offdiag.n_elem,2);
  dat2.col(0)=abs(vectorise(logcon(offdiag)));
  dat2.col(1)=vectorise(kinship(offdiag));
  int n2 = dat2.n_rows, k2 = dat2.n_cols;
  
  vec coef2 = solve(dat2, m_undir(offdiag)); 
  vec resid2 = m_undir(offdiag) - dat2*coef2; 
  
  double sig22 = as_scalar(trans(resid2)*resid2/(n2-k2));
  
  
  mat equated_rest=btransfers-logcon-theta(0)-theta(1)*kinship;
  vec sqyared_rest=pow(equated_rest(find(logcon>0)),2);
  
  
  vec sqyared_rest3=pow(equated_rest(offdiag),2);
  
  
  //double sqresidual_proxy;
  
//if (sqyared_rest.n_elem>0) {
double  sqresidual_proxy=mean(sqyared_rest);
double  sqresidual_proxy3=mean(sqyared_rest3);
//} else {
    //sqresidual_proxy=999;
//  }
  vec ret = {density,                   fb2,
             ib,                        sa,
             ra,                        pathlenghts,
             c1(0),                     c2(0),
             sig2,
             sig22,c3,sqresidual_proxy3};
  return(ret.replace(datum::nan,0));
}
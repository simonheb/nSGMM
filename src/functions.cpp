#include <RcppArmadillo.h>
#include <iostream>
//#include "gperftools/profiler.h"
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
mat timers=zeros(100,2);
int current=0;
auto start=std::chrono::system_clock::now();
void toc() {
  //  std::cout<<"ending" << current <<endl;
  timers(current,1)+=std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::system_clock::now() - start).count();
  current = 0;
  timers.col(0)=linspace(0,99);
}
void tic(int counter){ 
  toc();
  current = counter;
  //std::cout<<"starting" << current <<endl;
  start=std::chrono::system_clock::now();
}/*
 // [[Rcpp::export]]
 SEXP start_profiler(SEXP str) {
 ProfilerStart(Rcpp::as<const char*>(str));
 return R_NilValue;
 }
 
 // [[Rcpp::export]]
 SEXP stop_profiler() {
 ProfilerStop();
 return R_NilValue;
 }*/
// [[Rcpp::export]]
double  utlity_cppvec(const vec& consumption, const rowvec& altruism) {
  return as_scalar(altruism*log(consumption));
}

// [[Rcpp::export]]
double mynegutility_cpp(const vec& mytransfers, int i, mat transfers,const mat& altruism,const vec& income) {
  transfers.row(i-1) = trans(mytransfers);
  vec consumption = income + trans(sum(transfers,0))-sum(transfers,1);
  return -utlity_cppvec(consumption,altruism.row(i-1));
}

// [[Rcpp::export]]
double mynegutility_cpp_consumption(const vec& mytransfers, int i, mat transfers,const mat& altruism,vec consumption) {
  consumption -=  transfers.row(i-1).t() -  mytransfers;
  consumption(i-1) += accu(transfers.row(i-1)) - accu(mytransfers);
  return -utlity_cppvec(consumption,altruism.row(i-1));
}


uvec  rank_people(vec x,int i) {
  x(i-1)=-datum::inf;
  uvec indices=sort_index(x, "decrease");
  return indices.head(x.n_elem-1);
}


/*doesnt do capacities
// [[Rcpp::export]]
vec BBP_get_BR_analytically_cpp(uword i, mat transfers, const vec& income, const mat& altruism) {
  vec mytransfers = zeros(income.n_elem);
  vec oldtransfer = trans(transfers.row(i-1));
  transfers.row(i-1).zeros();
  
  vec consumption = income + trans(sum(transfers,0))-sum(transfers,1);
  vec marginal_util=altruism.col(i-1)/clamp(consumption,0.00000001,datum::inf);
  uvec corder=rank_people(marginal_util,i);
/////////////////////////////can i not be much faster if only take the index of the next best person evey round?
  vec suggestedtransfers;
  uvec me = {i-1};
  for(uword giveto=1; giveto<corder.n_elem; ++giveto) { 
    uvec whotogiveto = corder.head(giveto);
    double C=accu(consumption.elem(whotogiveto));
    double A=accu(altruism.submat(me,whotogiveto));
    double T=(C-A*consumption(i-1))/(-1-A);
    
    suggestedtransfers = altruism.submat(whotogiveto,me)*(consumption(i-1)-T)-consumption.elem(whotogiveto);
    if (suggestedtransfers.min()<0) {
      break;
    }
    mytransfers.elem(whotogiveto) = suggestedtransfers;
  }
  
  return mytransfers;
}*/
uword max_ind_simon(const vec & x) {
  double curr = -datum::inf;
  uword ret;
  for(uword i=0;i<x.n_elem;i++) {
    if (x(i)>curr) {
      ret=i;
      curr=x(i);
    }
  }
  return(ret);
}
bool BBP_update_BR_analytically_cpp_fast5(uword i, mat & transfers, vec & net_transfers_in, vec & net_transfers_out, const vec& income, const mat& altruism) {
  
  //set everything as if i was not giving anything
  vec mytransfers = zeros(income.n_elem);
  net_transfers_in -= trans(transfers.row(i-1));
  net_transfers_out(i-1)=0;
  
  vec consumption = income + net_transfers_in-net_transfers_out;
  vec marginal_util=altruism.col(i-1)/clamp(consumption,0.000000001,datum::inf); //the clamp is required because otherwise negative values or c (because of helpers dropping out) would imply negative marginal utlity
  
  vec suggestedtransfers;
  double Told=0;
  uvec me = {i-1};
  marginal_util(i-1)=-datum::inf;
  
  uvec whotogiveto = zeros<uvec>(income.n_elem);
  
  double C=0;
  double A=0;//accu(altruism.submat(me,whotogiveto.head(giveto)));
  
  for(uword giveto=1; giveto<income.n_elem; ++giveto) {  //if people have many transactions, it might make sense to not go through the first cases, this would however require a full ranking first
    whotogiveto(giveto-1) = max_ind_simon(marginal_util);
    
    marginal_util(whotogiveto(giveto-1))=-datum::inf;
    C+=consumption(whotogiveto(giveto-1)); // double C=accu(consumption.elem(whotogiveto.head(giveto)));
    A+=altruism(i-1,whotogiveto(giveto-1));// double A= accu(altruism.submat(me,whotogiveto.head(giveto)));
    double T=(C-A*consumption(i-1))/(-1-A);
    suggestedtransfers = altruism.submat(whotogiveto.head(giveto),me)*(consumption(i-1)-T)-consumption.elem(whotogiveto.head(giveto));
    if (suggestedtransfers.min()<0) { //possibly only need to test the last?
      break;
    }
    
    mytransfers.elem(whotogiveto.head(giveto)) = suggestedtransfers;
    Told=T;
  }
  bool updated=(!approx_equal(trans(mytransfers),transfers.row(i-1),"absdiff",0.00001));
  transfers.row(i-1)=trans(mytransfers);
  net_transfers_in += mytransfers; 
  net_transfers_out(i-1)=Told;
  return(updated);
}

bool BBP_update_BR_analytically_cpp_fast7_smarter(uword i, mat & transfers, vec & net_transfers_in, vec & net_transfers_out, const vec& income, const mat& altruism, const mat& capacities) {
  
  //set everything as if i was not giving anything
  vec mytransfers = zeros(income.n_elem);
  net_transfers_in -= trans(transfers.row(i-1));
  net_transfers_out(i-1)=0;
  
  vec consumption = income + net_transfers_in-net_transfers_out;
  vec marginal_util=altruism.col(i-1)/clamp(consumption,0.000000001,datum::inf); //the clamp is required because otherwise negative values or c (because of helpers dropping out) would imply negative marginal utlity
  
  vec suggestedinteriorsolutiontranstransfers;
  uvec me = {i-1}; //,
  marginal_util(i-1)=-datum::inf;
  
  uvec whotogiveto = zeros<uvec>(income.n_elem);
  uvec whotomaxout = zeros<uvec>(income.n_elem);
  uword n_maxout; //like giveto counts how many to give anything to, n_maxout counts how many to give max to
  
  double C=0;
  double A=0;//accu(altruism.submat(me,whotogiveto.head(giveto)));
  double T;

    
  for(uword giveto=1; giveto<income.n_elem; ++giveto) {  //if people have many transactions, it might make sense to not go through the first cases, this would however require a full ranking first
    whotogiveto(giveto-1) = max_ind_simon(marginal_util);
    n_maxout = 0;
    bool someoneisnotatcorner = true;
    bool therestisnotoptized = true;
    marginal_util(whotogiveto(giveto-1))=-datum::inf;
    double isconsumptionnetofmaxtransfers=consumption[i-1];
    C+=consumption(whotogiveto(giveto-1)); // double C=accu(consumption.elem(whotogiveto.head(giveto)));
    A+=altruism(i-1,whotogiveto(giveto-1));// double A= accu(altruism.submat(me,whotogiveto.head(giveto)));
    uvec whotooptimize=whotogiveto.head(giveto);
    double CC = C; //the one we use within this loop
    double AA = A; //-"- here we take out some obs every round
    while (someoneisnotatcorner&therestisnotoptized) {
      T=(CC-AA*(isconsumptionnetofmaxtransfers))/(-1-AA);
      suggestedinteriorsolutiontranstransfers = trans(altruism.submat(me,whotooptimize))*(isconsumptionnetofmaxtransfers-T)-consumption.elem(whotooptimize);
      therestisnotoptized=false;
      uvec excessive=find(suggestedinteriorsolutiontranstransfers>capacities(me,whotooptimize).t()); //who would receive too much?
      uvec newmaxout=whotooptimize(excessive);//their ids will be saved
      if (newmaxout.n_elem>0) { 
        therestisnotoptized=true;
        whotomaxout.subvec(n_maxout,n_maxout+newmaxout.n_elem-1)=newmaxout;
        n_maxout+=newmaxout.n_elem;
        for(uword persontomaxout: newmaxout) {
          //drop from those to optizie over
          whotooptimize.shed_row(conv_to<uword>::from(find(whotooptimize==persontomaxout)));
          
          //take out of CC and AA
          CC -= consumption(persontomaxout);
          AA -= altruism(i-1,persontomaxout);
          isconsumptionnetofmaxtransfers -= capacities(i-1,persontomaxout);
        }
        suggestedinteriorsolutiontranstransfers.reset();  //in case i move on now, (happens if all are maxed)
      }
      if (whotooptimize.n_elem==0) someoneisnotatcorner=false;
      
    }
    if (suggestedinteriorsolutiontranstransfers.n_elem>0){
      if (suggestedinteriorsolutiontranstransfers(suggestedinteriorsolutiontranstransfers.n_elem-1)<0) { //possibly only need to test the last? but that isn't any faster :/
        break;
      }
    } 
    mytransfers.elem(whotooptimize) = suggestedinteriorsolutiontranstransfers;
    mytransfers.elem(whotomaxout.head(n_maxout)) = capacities(me,whotomaxout.head(n_maxout));
  }
  bool updated=(!approx_equal(trans(mytransfers),transfers.row(i-1),"absdiff",0.00001));
  transfers.row(i-1)=trans(mytransfers);
  net_transfers_in += mytransfers; 
  net_transfers_out(i-1)=accu(mytransfers);
  return(updated);
}

// [[Rcpp::export]]
mat get_BBP_BR_analytically_cpp_inline(uword i, mat transfers, const vec& income, const mat& altruism, const mat& capacities) {
  std::cout<< 1<<endl;
  vec transfers_in=trans(sum(transfers,0));
  vec transfers_out=sum(transfers,1);
  BBP_update_BR_analytically_cpp_fast7_smarter(i,transfers,transfers_in,transfers_out,income,altruism,capacities);
  return(transfers.row(i-1));
}

// [[Rcpp::export]]
bool BBP_in_equilibrium_YaT_cpp(const mat& transfers, const vec& income, const mat& altruism,const mat& capacities, double tolerance) {
  vec consumption = income + trans(sum(transfers,0))-sum(transfers,1);
  mat con = (1/consumption)*trans(consumption);
  umat eqcondition = ((abs(con-altruism)<tolerance) % (transfers>0)) + ((transfers==0) % (con>=altruism));
  eqcondition.diag().ones();

  if (all(all(eqcondition))) {
    return true;
  }else{
    return false;
  }
}


// [[Rcpp::export]]
mat meannanrm(mat x, int dim) {
  mat x0=x;
  x0.elem(find_nonfinite(x)).zeros();
  mat x1=x0;
  x1.elem(find_finite(x)).ones();
  return(sum(x0,dim)/sum(x1,dim));
}

// [[Rcpp::export]]
mat consumption_weights_cpp(const mat& alphas,  const mat transfers, const mat& t_conmat) {
  int n=transfers.n_rows;
  mat transferdirections=transfers-trans(transfers);
  transferdirections.replace(0,datum::nan);
  mat care=exp(transferdirections%log(alphas)); //pow(alphas, transferdirections);
  care.replace(datum::nan,0);
  mat care0=care;
  
  mat cw=zeros(n,n);
  cw.replace(0,datum::nan);
  
  for(int zz=0; zz<n; ++zz) {
    uvec l = find_nonfinite(cw);
    if (max(care(l))==0) {break;} //there's nothing to get anymore
    cw(l)=care(l);
    cw.replace(0,datum::nan);
    care=care * care0;
    care.diag().zeros();
  }
  cw.diag().ones(); //CW tells us about how consumptions must relate to each other to fit the equality constraint for linked hhs
  cw(find(t_conmat==0)).fill(datum::nan);
  vec consumption_fractions=trans(meannanrm(cw.each_col() / meannanrm(cw,1),0));
  return consumption_fractions;
}


// [[Rcpp::export]]
vec BBP_c_from_atY_cpp(const mat& alphas,const mat& transferstructure,const vec& incomes,const mat& t_components_matrix, const uvec& t_components_csize,const mat& t_conmat) {
  
  //get component-wise average incomes
  vec component_incomes = trans(incomes.t()*t_components_matrix)  / t_components_csize;
  //get within-component distribution from alphas (only works via tree structure)
  vec cw=consumption_weights_cpp(alphas,transferstructure,t_conmat);
  mat full_incomes_by_comp = t_components_matrix*component_incomes;
  
  return ( full_incomes_by_comp % cw);
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
  //loop over is
  //int *array_of_pointers[transferstructure.n_rows]; 
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
mat BBP_T_from_tYc_cpp(mat transferstructure, const vec& incomes, const vec& consumptions, mat Tr, uword depth=0) {
  //probably this can be significantly sped up by refusing to wokr on cases where transfers are inconsistent with consumption
  //find nodes that have only one t
  
  if (depth==0) {
    Tr = zeros(transferstructure.n_rows,transferstructure.n_rows);
  }
  //receivers
  vec indegree=trans(sum(transferstructure,0));
  vec outdegree=sum(transferstructure,1);
  vec intermediateincomes=incomes+trans(sum(Tr,0))-sum(Tr,1);
  uvec canTupdate=find((indegree!=1) + (outdegree!=0));
  vec netflows_to_receivers=consumptions-intermediateincomes; //otherwise this counts double
  
  if(any((indegree>0) % (outdegree==0) % (netflows_to_receivers<0))) {
    //std::cout << "not an eq"<< endl;
    return Tr;
  }
  
  netflows_to_receivers(canTupdate).zeros();
  
  mat newTr=transferstructure.each_row()%trans(netflows_to_receivers);
  
  Tr(find(newTr>0))=newTr(find(newTr>0));
  transferstructure(find(Tr>0)).zeros();
  
  //givers
  
  indegree=trans(sum(transferstructure,0));
  outdegree=sum(transferstructure,1);
  intermediateincomes=incomes+trans(sum(Tr,0))-sum(Tr,1);
  canTupdate=find((outdegree!=1) + (indegree!=0));
  
  vec netflows_from_givers=consumptions-intermediateincomes; //otherwise this counts double
  
  
  if(any((outdegree>0) % (indegree==0) % (netflows_from_givers>0))) {
    //std::cout << "Cnot an eqqwe"<< endl;
    return Tr;
  }
  
  netflows_from_givers(canTupdate).zeros();
  newTr=transferstructure.each_col()%(-netflows_from_givers);
  Tr(find(newTr>0))=newTr(find(newTr>0));
  transferstructure(find(Tr>0)).zeros();
  
  if(accu(transferstructure)==0) { //we're done
    return Tr;
  } else if(depth>transferstructure.n_rows) {
    //  std::cout << "cant help, this happens if we cant solve ti" << endl;
    return Tr;
  } else {
    return BBP_T_from_tYc_cpp(transferstructure,incomes,consumptions,Tr,depth=depth+1);
  }
}
// [[Rcpp::export]]
mat BBP_T_from_atY_plain_cpp(const mat& alphas, const mat& transferstructure, const vec& incomes) {
  
  if (any(transferstructure(find(transferstructure>0))<1)) {
    std::cout << "this shouldnt be non-binary" << transferstructure << endl;;
  }
  
  uvec t_components_csize;
  uvec t_components_membership;
  mat t_components_matrix;
  mat t_conmat;
  Ccomponents(transferstructure, t_components_csize, t_components_membership, t_components_matrix, t_conmat);
  
  
  
  vec consumptions = BBP_c_from_atY_cpp(alphas,transferstructure,incomes,t_components_matrix,t_components_csize,t_conmat);
  
  return BBP_T_from_tYc_cpp(transferstructure,incomes,consumptions,zeros(alphas.n_rows,alphas.n_cols));
}

/* this does not take car of capacities
// [[Rcpp::export]]
mat equilibrate_cpp(const mat& altruism, const vec& income,int modmode=5) {
  mat transfers = zeros(income.n_elem,income.n_elem);
  int updates;
  int r;
  for (r=0;r<20000;r++) {
    if (r%20==modmode) {    //lets try if the structure is already final, in which case we can find the ne algebraically
      
      mat transferstructure=transfers;
      transferstructure.replace(0,datum::nan);
      transferstructure(find_finite(transferstructure)).ones();
      transferstructure.replace(datum::nan,0);
      mat TC = BBP_T_from_atY_plain_cpp(altruism,transferstructure,income);
      //std::cout << all(all(TCa==TC)) << endl;
      
      bool eqq = BBP_in_equilibrium_YaT_cpp(TC, income,altruism,capacity,0.00001);
      if(eqq) {
        transfers = TC;
      }
      
    }
    updates=0;
    //std::cout << "C" << endl;
    uvec updating = linspace<uvec>(1, income.n_elem,income.n_elem);
    uword updatingmax = income.n_elem;
    for (int rep=0;rep<1;rep++){
      uword next_updateingmax = 0;
      //std::cout <<endl << "("<<rep<<"): ";
      for (uword u=0; u<updatingmax; u++) { //for each player
        uword i=updating(u);        //find the best response
        //std::cout << endl<<i;
        vec BR=BBP_get_BR_analytically_cpp(i,transfers,income,altruism,capacity);
        //std::cout << trans(BR);
        //update transfers
        if (!approx_equal(transfers.row(i-1),trans(BR),"absdiff",0.00001)) {
          //std::cout <<  "updated, ";
          updates=updates+1;
          updating(next_updateingmax++)=i;
        }
        transfers.row(i-1) = trans(BR);
      }
      if (next_updateingmax==0) break;
      updatingmax = next_updateingmax;
    }
//    if ((r%20==0)|(updates==0)) {std::cout << "c:Round" << r <<"... (" <<updates <<" nodes updated their transactions)" << endl;}
    if (updates==0) {break;}
  }
  //std::cout <<endl << r<<endl;
  transfers=max(transfers-trans(transfers),zeros(income.n_elem,income.n_elem));
  if (updates>0) {
    std::cout<< "c:Best responses did not converge to a NE, probably you need to increase the rounds" << endl;
    return transfers;
  }
  //  else {std::cout << "F:Stopped after "<< r <<" rounds. Found a/the nash equilibium" << updates << endl;}
  return transfers ;
}
*/
 
 
 // [[Rcpp::export]]
 mat equilibrate_cpp_fast5(const mat& altruism, const vec& income,int modmode=5) {
   mat transfers = zeros(income.n_elem,income.n_elem);
   vec net_transfers_in  = zeros(income.n_elem);
   vec net_transfers_out = zeros(income.n_elem);
   int updates; 
   int r;
   for (r=0;r<30000;r++) {
     updates=0;
     uvec updating = linspace<uvec>(1, income.n_elem,income.n_elem);
     uword updatingmax = income.n_elem;
     for (int rep=0;rep<20;rep++){
       uword next_updateingmax = 0;
       //std::cout <<endl << "("<<rep<<"): ";
       for (uword u=0; u<updatingmax; u++) { //for each player
         uword i=updating(u);        //find the best response
         //std::cout << i;
         bool updated=BBP_update_BR_analytically_cpp_fast5(i,transfers,net_transfers_in,net_transfers_out,income,altruism);
         //update transfers
         if (updated) {//!approx_equal(transfers.row(i-1),oldt,"absdiff",0.00001)) {
           updates=updates+1;
           updating(next_updateingmax++)=i;
         }
       }
       if (next_updateingmax==0) break;
       updatingmax = next_updateingmax;
     }
     //if ((r%20==0)|(updates==0)|(1==1)) {std::cout << "c:Round" << r <<"... (" <<updates <<" nodes updated their transactions)" << endl;}
     if (updates==0) {break;}
   }
   //  std::cout << timers.head_rows(3) << endl; 
   transfers=max(transfers-trans(transfers),zeros(income.n_elem,income.n_elem));
   if (updates>0) {
     std::cout<< "c:Best responses did not converge to a NE, probably you need to increase the rounds" << endl;
     return transfers;
   } 
   //std::cout << "r="<<r << endl;
   //  else {std::cout << "F:Stopped after "<< r <<" rounds. Found a/the nash equilibium" << updates << endl;}
   return transfers ;
 }


/*
// [[Rcpp::export]]
mat equilibrate_cpp_fast5_smarter(const mat& altruism, const vec& income, const mat& capacity,int modmode=5) {
  mat transfers = zeros(income.n_elem,income.n_elem);
  vec net_transfers_in  = zeros(income.n_elem);
  vec net_transfers_out = zeros(income.n_elem);
  int updates; 
  int r;
  for (r=0;r<30000;r++) {
    updates=0;
    uvec updating = linspace<uvec>(1, income.n_elem,income.n_elem);
    uword updatingmax = income.n_elem;
    for (int rep=0;rep<20;rep++){
      uword next_updateingmax = 0;
      //std::cout <<endl << "("<<rep<<"): ";
      for (uword u=0; u<updatingmax; u++) { //for each player
        uword i=updating(u);        //find the best response
        //std::cout << i;
        bool updated=BBP_update_BR_analytically_cpp_fast5_smarter(i,transfers,net_transfers_in,net_transfers_out,income,altruism,capacity);
        //update transfers
        if (updated) {//!approx_equal(transfers.row(i-1),oldt,"absdiff",0.00001)) {
          updates=updates+1;
          updating(next_updateingmax++)=i;
        }
      }
      if (next_updateingmax==0) break;
      updatingmax = next_updateingmax;
    }
    //if ((r%20==0)|(updates==0)|(1==1)) {std::cout << "c:Round" << r <<"... (" <<updates <<" nodes updated their transactions)" << endl;}
    if (updates==0) {break;}
  }
  //  std::cout << timers.head_rows(3) << endl; 
  transfers=max(transfers-trans(transfers),zeros(income.n_elem,income.n_elem));
  if (updates>0) {
    std::cout<< "c:Best responses did not converge to a NE, probably you need to increase the rounds" << endl;
    return transfers;
  } 
  //std::cout << "r="<<r << endl;
  //  else {std::cout << "F:Stopped after "<< r <<" rounds. Found a/the nash equilibium" << updates << endl;}
  return transfers ;
}

*/

// [[Rcpp::export]]
mat equilibrate_cpp_fast7_smarter(const mat& altruism, const vec& income, const mat& capacity,int modmode=5) {
  mat transfers = zeros(income.n_elem,income.n_elem);
  vec net_transfers_in  = zeros(income.n_elem);
  vec net_transfers_out = zeros(income.n_elem);
  int updates; 
  int r;
  for (r=0;r<30000;r++) {
    updates=0;
    uvec updating = linspace<uvec>(1, income.n_elem,income.n_elem);
    uword updatingmax = income.n_elem;
    for (int rep=0;rep<20;rep++){
      uword next_updateingmax = 0;
      //std::cout <<endl << "("<<rep<<"): ";
      for (uword u=0; u<updatingmax; u++) { //for each player
        uword i=updating(u);        //find the best response
        //std::cout << i;
        bool updated=BBP_update_BR_analytically_cpp_fast7_smarter(i,transfers,net_transfers_in,net_transfers_out,income,altruism,capacity);
        //update transfers
        if (updated) {//!approx_equal(transfers.row(i-1),oldt,"absdiff",0.00001)) {
          updates=updates+1;
          updating(next_updateingmax++)=i;
        }
      }
      if (next_updateingmax==0) break;
      updatingmax = next_updateingmax;
    }
    //if ((r%20==0)|(updates==0)|(1==1)) {std::cout << "c:Round" << r <<"... (" <<updates <<" nodes updated their transactions)" << endl;}
    if (updates==0) {break;}
  }
  //  std::cout << timers.head_rows(3) << endl; 
  transfers=max(transfers-trans(transfers),zeros(income.n_elem,income.n_elem));
  if (updates>0) {
    std::cout<< "c:Best responses did not converge to a NE, probably you need to increase the rounds" << endl;
    return transfers;
  } 
  //std::cout << "r="<<r << endl;
  //  else {std::cout << "F:Stopped after "<< r <<" rounds. Found a/the nash equilibium" << updates << endl;}
  return transfers ;
}


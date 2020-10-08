#include <RcppArmadillo.h>
#include <iostream>
#include <tictoc.h>
#include <network_functions.h>
using namespace arma;



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



int max_ind_simon(const vec & x) {
  double curr = -datum::inf;
  uword ret=-1;
  for(uword i=0;i<x.n_elem;i++) {
    if (x(i)>curr) {
      ret=i;
      curr=x(i);
    }
  }
  return(ret);
}

bool BBP_update_BR_analytically_cpp_fast10(uword i, mat & transfers, vec & net_transfers_in, vec & net_transfers_out, const vec& income, const mat& altruism, const mat& capacities) {
  //set everything as if i was not giving anything
  vec mytransfers = zeros(income.n_elem);
  net_transfers_in -= trans(transfers.row(i-1));
  net_transfers_out(i-1)=0;

  vec consumption = income + net_transfers_in-net_transfers_out;
  vec marginal_util=altruism.row(i-1).t()/clamp(consumption,0.000000001,datum::inf); //the clamp is required because otherwise negative values or c (because of helpers dropping out) would imply negative marginal utlity

  vec suggestedinteriorsolutiontranstransfers;
  uvec me = {i-1}; //,
  marginal_util(i-1)=-datum::inf;

  uvec whotogiveto = zeros<uvec>(income.n_elem);
  uvec whotomaxout = zeros<uvec>(income.n_elem);
  uword n_maxout; //like giveto counts how many to give anything to, n_maxout counts how many to give max to
  
  double C=0;
  double A=0;//accu(altruism.submat(me,whotogiveto.head(giveto)));

  int nextrecipient = max_ind_simon(marginal_util);
  marginal_util(nextrecipient)=-datum::inf;
  double nextC=consumption(nextrecipient);
  double nextA=altruism(i-1,nextrecipient);
  for(uword giveto=1; giveto<income.n_elem; giveto++) {  //if people have many transactions, it might make sense to not go through the first cases, this would however require a full ranking first
    n_maxout = 0;
    whotogiveto(giveto-1) = nextrecipient;
    //pick the one to use next round  
    nextrecipient = max_ind_simon(marginal_util);
    double isconsumptionnetofmaxtransfers=consumption[i-1];
    if (nextrecipient!=-1) {
      marginal_util(nextrecipient)=-datum::inf;

      C=nextC;
      nextC+=consumption(nextrecipient); // double C=accu(consumption.elem(whotogiveto.head(giveto)));
      A=nextA;
      nextA+=altruism(i-1,nextrecipient);// double A= accu(altruism.submat(me,whotogiveto.head(giveto)));
      
      //if adding mone more recipient still does not lead to negative transfers
      if ( (isconsumptionnetofmaxtransfers+nextC)/consumption(nextrecipient) > (1+nextA)/altruism(i-1,nextrecipient)) {continue;}
    }
    uvec whotooptimize=whotogiveto.head(giveto);
    double CC = C; //the one we use within this loop
    double AA = A; //-"- here we take out some obs every round
    
    bool someoneisnotatcorner = true;
    bool therestisnotoptized = true;
    while (someoneisnotatcorner&therestisnotoptized) {
      suggestedinteriorsolutiontranstransfers = trans(altruism.submat(me,whotooptimize))*(isconsumptionnetofmaxtransfers-(CC-AA*(isconsumptionnetofmaxtransfers))/(-1-AA))-consumption.elem(whotooptimize);
      therestisnotoptized=false;
      uvec excessive=find(suggestedinteriorsolutiontranstransfers>capacities(me,whotooptimize).t()); //who would receive too much?
      uvec newmaxout=whotooptimize(excessive);//their ids will be saved
      if (newmaxout.n_elem>0) {
        whotomaxout.subvec(n_maxout,n_maxout+newmaxout.n_elem-1)=newmaxout;
        n_maxout+=newmaxout.n_elem;
        if (newmaxout.n_elem==whotooptimize.n_elem) {
          someoneisnotatcorner=false;
          suggestedinteriorsolutiontranstransfers.reset();  //in case i move on now, (happens if all are maxed)
          whotooptimize.reset();
        } else {
          therestisnotoptized=true;
          for(uword persontomaxout: newmaxout) {
            //drop from those to optizie over
            whotooptimize.shed_row(conv_to<uword>::from(find(whotooptimize==persontomaxout)));
            
            //take out of CC and AA
            CC -= consumption(persontomaxout);
            AA -= altruism(i-1,persontomaxout);
            isconsumptionnetofmaxtransfers -= capacities(i-1,persontomaxout);
          }
        }

      }
      
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
  BBP_update_BR_analytically_cpp_fast10(i,transfers,transfers_in,transfers_out,income,altruism,capacities);
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


// [[Rcpp::export]]
mat equilibrate_cpp_fast8_debug(const mat& altruism, const vec& income, const mat& capacity,bool verbose, int& r, bool& NE, int maxrounds=500) {
  mat transfers = zeros(income.n_elem,income.n_elem);
  vec net_transfers_in  = zeros(income.n_elem);
  vec net_transfers_out = zeros(income.n_elem);
  int updates; 
//std::cout << "(" << altruism(1,2) <<  ")";
  for (r=0;r<maxrounds;r++) {
    updates=0;
    uvec updating = linspace<uvec>(1, income.n_elem,income.n_elem);
    uword updatingmax = income.n_elem;
    for (int rep=0;rep<20;rep++){
      uword next_updateingmax = 0;
      //std::cout <<endl << "("<<rep<<"): ";
      for (uword u=0; u<updatingmax; u++) { //for each player
        uword i=updating(u);        //find the best response
        //std::cout << i;
        bool updated=BBP_update_BR_analytically_cpp_fast10(i,transfers,net_transfers_in,net_transfers_out,income,altruism,capacity);
        //update transfers
        if (updated) {//!approx_equal(transfers.row(i-1),oldt,"absdiff",0.00001)) {
          updates=updates+1;
          updating(next_updateingmax++)=i;
        }
      }
      if (next_updateingmax==0) break;
      updatingmax = next_updateingmax;
    }
//  if (verbose & ((r%20==0)|(updates==0)|(1==1))) {std::cout << "c:Round" << r <<"... (" <<updates <<" nodes updated their transactions)" << endl;}
    if (updates==0) {break;}
  }
  
  //std::cout << "r:" <<r;
  transfers=max(transfers-trans(transfers),zeros(income.n_elem,income.n_elem));
  //if (verbose) std::cout << endl << r<<endl;
  //if (verbose & r>200) std::cout << "********";
  NE=true;
  if (updates>0) {
    NE=false;
    //std::cout<< "c:Best responses did not converge to a NE (yet)" << endl;
    return transfers;
  } 
  return transfers ;
}


// [[Rcpp::export]]
mat equilibrate_cpp_fast8_smarter(const mat& altruism, const vec& income, const mat& capacity, int maxrounds=500) {
  int r;
  bool NE;
  return(equilibrate_cpp_fast8_debug(altruism, income, capacity,false,r,NE));
}
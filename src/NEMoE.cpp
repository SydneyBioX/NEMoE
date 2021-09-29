#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(cpp11)]]

/* -------------------------------------------
 ---------- Untils Function ------------
 ------------------------------------------- */

double softThresh(double x, double lambda) {

  double y;
  if(x > lambda){
    y = x - lambda;
  }else if(x < (-lambda)){
    y = x + lambda;
  }else{
    y = 0;
  }
  return y;
}

double clipping(double x, double x_min, double x_max){
  if(x < x_min){
    x = x_min;
  }else if(x > x_max){
    x = x_max;
  }
  return x;
}

NumericVector clippingV(Rcpp::NumericVector x, double x_min, double x_max){

  int n = x.length();
  NumericVector y = x;
  for(int i = 0; i<n; i++){
    y[i] = clipping(x[i], x_min, x_max);
  }
  return y;
}

Eigen::VectorXd clippingVX(Eigen::VectorXd x, double x_min, double x_max){

  int n = x.size();
  for(int i = 0; i<n; i++){
    x(i) = clipping(x(i), x_min, x_max);
  }
  return x;
}

Eigen::MatrixXd clippingMX(Eigen::MatrixXd x, double x_min, double x_max){

  int p = x.cols();
  for(int i = 0; i<p; i++){
    x.col(i) = clippingVX(x.col(i), x_min, x_max);
  }
  return x;
}

NumericVector checking(NumericVector x, double na_map = 0,
                       double inf_map = 1e3){

  LogicalVector x_na, x_inf;
  int n = x.length();
  x_na = is_nan(x);
  x_inf = is_infinite(x);
  for(int i = 0; i <n; i++){
    if(x_na(i)){
      x(i) = na_map;
    }else if(x_inf(i)){
      x(i) = inf_map;
    }
  }
  return x;
}

NumericVector softmax(NumericVector x){
  NumericVector res;
  res = exp(x)/sum(exp(x));
  return res;
}

Eigen::MatrixXd logit(Eigen::MatrixXd X, Eigen::VectorXd beta){

  Eigen::VectorXd odd = X * beta;
  NumericVector odd_1 = wrap(odd);
  NumericVector prob;
  prob = 1/(1 + exp(-odd_1));
  Eigen::MatrixXd prob0;
  prob0 = as<Eigen::Map<Eigen::VectorXd>>(prob);
  return prob0;
}

Eigen::MatrixXd Exccol(Eigen::MatrixXd A, int r){

  int n = A.rows(), p = A.cols();
  Eigen::MatrixXd B(n, p - 1);
  if(r == 0){
    B = A.rightCols(p - 1);
  }else if(r == (p - 1)){
    B = A.leftCols(p - 1);
  }else{
    B.leftCols(r) = A.leftCols(r);
    B.rightCols(p - r - 1) = A.rightCols(p - r - 1);
  }

  return B;
}

double pen_V(Eigen::VectorXd beta, NumericVector lambda, double alpha){

  int p = lambda.size();
  NumericVector beta1 = wrap(beta.tail(p));
  double pen = 0.5*(1- alpha)*sum(lambda*pow(beta1, 2)) +
    alpha*sum(lambda*abs(beta1));
  return pen;
}

double pen_M(Eigen::MatrixXd Beta, NumericVector lambda, double alpha){
  int K = Beta.cols();
  double pen = 0;
  for(int k =0;k<K;k++){
    pen += pen_V(Beta.col(k), lambda, alpha);
  }
  return pen;
}

IntegerVector hard_idx(NumericMatrix X){

  int n = X.rows(), p = X.cols();
  IntegerVector Idx(n);

  for(int i = 0; i<n; i++){
    if(p == 1){
      Idx[i] = round(X(i,0));
    }else{
      Idx[i] = which_max(X.row(i));
    }
  }
  return Idx;
}

IntegerMatrix hard_prob(NumericMatrix X){
  int n = X.rows(), p = X.cols();
  IntegerMatrix H(n,p);
  H.fill(0);
  IntegerVector Idx = hard_idx(X);
  for(int i=0; i<n; i++){
    H(i, Idx[i]) = 1;
  }
  return H;

}

IntegerVector pos1(IntegerVector x){
  int n = x.size(), k = 0;
  IntegerVector pos(sum(x));
  for(int i =0; i<n; i++){
    if(x[i] == 1){
      pos[k] = i;
      k = k + 1;
    }
  }
  return pos;
}

NumericMatrix RowSel(NumericMatrix X, IntegerVector idx){

  int p = X.cols(), n = idx.size();
  NumericMatrix Y(n,p);
  for (int i =0; i<n; i++){
    Y.row(i) = X.row(idx[i]);
  }

  return Y;
}

Eigen::MatrixXd RowSelX(Eigen::MatrixXd X, IntegerVector idx){

  NumericMatrix X1, Y1;
  X1 = wrap(X);
  Y1 = RowSel(X1, idx);
  Eigen::MatrixXd Y = as<Eigen::Map<Eigen::MatrixXd>>(Y1);

  return Y;
}

NumericMatrix ColSel(NumericMatrix X, IntegerVector idx){

  int n = X.rows(), p = idx.size();
  NumericMatrix Y(n,p);
  for (int i =0; i<p; i++){
    Y.column(i) = X.column(idx[i]);
  }

  return Y;
}

Eigen::MatrixXd ColSelX(Eigen::MatrixXd X, IntegerVector idx){

  NumericMatrix X1, Y1;
  X1 = wrap(X);
  Y1 = ColSel(X1, idx);
  Eigen::MatrixXd Y = as<Eigen::Map<Eigen::MatrixXd>>(Y1);

  return Y;
}

NumericMatrix R_sum(NumericMatrix r_i, int L){

  int nL = r_i.rows(), K = r_i.cols(), n = nL/L;
  NumericMatrix r_i_s(n, K);
  IntegerVector slice = seq(0, L - 1)*n;
  for(int i = 0; i< n; i++){
    r_i_s.row(i) = colMeans(RowSel(r_i, slice + i));
  }
  return r_i_s;
}

//' Cpp function: Sample from a Dirichlet distribution
//' @param n Number of random vector to return.
//' @param a A positive vector of parameters.
//' @return A random matrix with n row from Dirichlet
//'  distribution with parameter a.
// [[Rcpp::export]]
NumericMatrix rDirichlet(int n, NumericVector a){
  int l = a.size();
  NumericMatrix x(n,l);
  NumericVector x_s(n);
  for(int i =0; i<l; i++){
    x.column(i) = Rcpp::rgamma(n, a[i]);
  }
  x_s = rowSums(x);
  for(int i = 0; i<l; i++){
    x.column(i) = x.column(i)/x_s;
  }
  return x;
}

//' Cpp function: Sample from Categorical distribution
//'
//' @param Prob A matrix with n row and p column. Each column represent a class.
//' The ith row is the probability distribution of ith sample.
//' @return A matrix of n row and p column. The ith row draw one
//' sample from categorical distribution with probability in the ith row in Prob.
// [[Rcpp::export]]
NumericMatrix rSample(NumericMatrix Prob){

  int n = Prob.rows(), K = Prob.cols();
  NumericMatrix X(n, K);
  for(int i =0;i < n; i++){
    IntegerVector X_i(n);
    NumericVector Prob_i = Prob.row(i);
    rmultinom(1, Prob_i.begin(), K, X_i.begin());
    X.row(i) = X_i;
  }
  return X;
}

double f_proj(NumericVector x, double alpha, double y){
  double f = sum(0.5*(1 - alpha)* pow(x - y, 2) + alpha*abs(x - y));
  return f;
}

double innerSearch(NumericVector beta1, double alpha){
  int idx;
  double y, beta1_med, beta1_mean, m, M;
  beta1_med = median(beta1);
  beta1_mean = mean(beta1);
  NumericVector mM = {beta1_med, beta1_mean};
  m = min(mM);
  M = max(mM);
  NumericVector f_eval(51), y_temp(51);

  if(alpha == 1){
    y = beta1_med;
  }else if(alpha == 0){
    y = beta1_mean;
  }else{

    for(int i = 0; i <= 50; i++){
      y_temp[i] = m + double(i)/50*(M - m);
      f_eval[i] = f_proj(beta1, alpha, y_temp[i]);
    }
    idx = which_min(f_eval);
    y = y_temp[idx];

  }
  return(y);

}

// [[Rcpp::export]]
Eigen::MatrixXd projMultiConstraint(Eigen::MatrixXd Beta0, double alpha = 1){

  Eigen::MatrixXd Beta = Beta0;
  NumericVector beta1;
  double c;
  int p = Beta0.rows();
  for(int k =0; k < p; k++){
    beta1 = wrap(Beta0.row(k));
    if(k == 0){
      beta1 = beta1 - mean(beta1);
    }else{
      c = innerSearch(beta1, alpha);
      beta1 = beta1 - c;
    }
    Beta.row(k) = as<Eigen::Map<Eigen::MatrixXd>>(beta1);
  }
  return Beta;

}

int calc_nz_row(Eigen::VectorXd W, double thresh = 1e-8){
  int p = W.size();
  int nz = 0;
  for(int i =1; i<p; i++){
    if(abs(W[i]) > thresh){
      nz ++;
    }
  }
  return nz;
}

/* -------------------------------------------
 ---------- Likelihood Function ------------
 ------------------------------------------- */

NumericVector LBinom0(NumericVector y, NumericVector prob, bool logout = false){

  int n = y.size();
  NumericVector LL(n);
  Eigen::ArrayXd y1(n), prob1(n);
  y1 = as<Eigen::Map<Eigen::ArrayXd>>(y);
  prob1 = as<Eigen::Map<Eigen::ArrayXd>>(prob);
  prob = clippingV(prob, 1e-9, 1 - 1e-9);
  if(logout){
    LL = y*log(prob) + (1 - y)*log(1 - prob);
  }else{
    LL = wrap((prob1.pow(y1))*((1 - prob1).pow(1 - y1)));
  }

  return LL;
}

NumericVector LBinom(Eigen::MatrixXd X, Eigen::VectorXd beta, NumericVector y,
                     NumericVector residual_k, bool logout = false){

  Eigen::VectorXd odd = X * beta;
  NumericVector odd_1 = wrap(odd);
  NumericVector prob;
  prob = 1/(1 + residual_k * exp(-odd_1));
  NumericVector LL;
  LL = LBinom0(y, prob, logout);
  return LL;
}

//' Cpp function: Calculate the estimated probability using softmax function
//' (multinomial regression).
//'
//' @param X a data matrix of input in multinomial regression.
//' @param B coefficients in multinomial regression.
//'
//' @return A matrix of multinomial probability with
//' p_i=exp(X %*% B_i)/sum(X %*% B)
//'
// [[Rcpp::export]]
Eigen::MatrixXd calcProb(Eigen::MatrixXd X, Eigen::MatrixXd B){

  int n = X.rows(), k = B.cols();
  Eigen::MatrixXd odd, prob(n, k);

  odd = X*B;

  for(int i = 0; i < n; i++){
    NumericVector odd_i = wrap(odd.row(i));
    NumericVector prob_i = softmax(odd_i);
    prob.row(i) = as<Eigen::Map<Eigen::VectorXd>>(prob_i);
  }
  return prob;
}

NumericVector LMulti0(Eigen::MatrixXd Y, Eigen::MatrixXd Prob,
                      bool logout = false){

  int n = Y.rows(), p = Y.cols();
  Eigen::ArrayXXd pow_res(n, p);
  Eigen::VectorXd LL1(n);
  NumericVector LL(n);
  pow_res = Prob.array().pow(Y.array());
  if(logout){
    pow_res = log(pow_res);
    LL1 = pow_res.matrix().rowwise().sum();
  }else{
    LL1 = pow_res.matrix().rowwise().prod();
  }
  LL = wrap(LL1);
  return LL;
}

NumericVector LMulti(Eigen::MatrixXd X, Eigen::MatrixXd B,
                       Eigen::MatrixXd Y, bool logout = false){

  Eigen::MatrixXd Prob = calcProb(X,B);

  NumericVector LL = LMulti0(Y, Prob, logout);

  return LL;
}


//' Cpp function: Predict probability of mixture distribution
//' @param X a data matrix of input in experts network.
//' @param Z a data matrix of input in gating network.
//' @param y a vector of response.
//' @param W parameters in experts network.
//' @param V parameters in gating network.
//' @return A matrix of predicted probability with
//' predProb = pi * (1/(1 + exp(-X*W))),
//' pi = exp(Z * V_i)/sum(Z * V)
// [[Rcpp::export]]
Eigen::MatrixXd predProb(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericVector y,
                          Eigen::MatrixXd W, Eigen::MatrixXd V){
  int n = X.rows(), p = X.cols(), K = V.cols(), q = Z.cols();
  Eigen::MatrixXd X1(n, p+1), Z1(n, q+1), ones = MatrixXd::Ones(n,1);
  X1 << ones, X; Z1 << ones, Z;
  Eigen::MatrixXd prob(n, K);
  Eigen::MatrixXd latent_prob = calcProb(Z1, V);
  Eigen::MatrixXd sub_prob(n,K);
  for(int k=0; k<K; k++){
    sub_prob.col(k) = logit(X1, W.col(k));
  }
  prob = (latent_prob.array() * sub_prob.array()).matrix();
  return prob;
}

NumericVector LL_obs0(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericVector y,
                      Eigen::MatrixXd W, Eigen::MatrixXd V,
                      bool logout = false){
  int n = X.rows(), K = V.cols();
  Eigen::VectorXd prob, onek = MatrixXd::Ones(K,1);
  NumericVector prob1(n), Likelihood(n);
  prob = predProb(X, Z, y, W, V) * onek;
  prob1 = wrap(prob);
  Likelihood = LBinom0(y, prob1);

  if(logout){
    Likelihood = clippingV(Likelihood, 1e-9, 1 - 1e-9);
    Likelihood = log(Likelihood);
  }
  return Likelihood;
}

NumericVector LL_obs(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericMatrix y,
                     NumericVector seg, Eigen::MatrixXd W, Eigen::MatrixXd V){

  int n = X.rows(), L = seg.size(), K = V.cols();

  NumericVector start_X(L), end_X(L), start_W(L), end_W(L), LL(L + 1);
  double LL_s = 0, LL_temp = 0;
  for(int l=0; l<L;l++){
    if(l == 0){
      start_X[l] = 0;
      start_W[l] = 0;
      end_X[l] = seg[l];
      end_W[l] = seg[l] + 1;
    }else{
      start_X[l] = end_X[l - 1];
      start_W[l] = end_W[l - 1];
      end_X[l]  = start_X[l] + seg[l];
      end_W[l]  = start_W[l] + seg[l] + 1;
    }
    LL_temp = sum(LL_obs0(X.block(0, start_X[l], n, seg[l]), Z, y.column(l),
                          W.block(start_W[l], 0, seg[l] + 1, K), V, true));

    LL[l] = LL_temp;
    LL_s += LL_temp;
  }
  LL[L] = LL_s;
  return LL;
}

NumericVector Pen_L(NumericVector seg, Eigen::MatrixXd W, Eigen::MatrixXd V,
                    NumericMatrix lambda1, NumericVector lambda2, int n,
                    NumericMatrix alpha1, double alpha2, bool multiply = false){

  int L = seg.size(), K = W.cols();
  NumericMatrix lambda11 = clone(lambda1);
  NumericVector lambda21 = clone(lambda2);
  if(multiply){
    lambda11 = lambda1 * n;
    lambda21 = lambda2 * n;
  }

  double pen2 = pen_M(V, lambda21, alpha2), pen1_s = 0;
  NumericVector start_X(L), end_X(L), start_W(L), end_W(L), pen1(L + 1);
  for(int l=0; l<L;l++){
    double pen1_temp = 0;
    if(l == 0){
      start_X[l] = 0;
      start_W[l] = 0;
      end_X[l] = seg[l];
      end_W[l] = seg[l] + 1;
    }else{
      start_X[l] = end_X[l - 1];
      start_W[l] = end_W[l - 1];
      end_X[l]  = start_X[l] + seg[l];
      end_W[l]  = start_W[l] + seg[l] + 1;
    }
    Eigen::MatrixXd W_temp = W.block(start_W[l], 0, seg[l] + 1, K);
    NumericMatrix lambda_temp = lambda11(seq(start_X[l],end_X[l] - 1),_);

    for(int k=0; k<K; k++){
      pen1_temp += pen_V(W_temp.col(k), lambda_temp.column(k), alpha1(l,k));
    }
    pen1[l] = pen1_temp + pen2;
    pen1_s += pen1_temp;
  }
  pen1[L] = pen1_s + pen2;

  return pen1;
}

NumericVector PLL_obs(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericMatrix y,
                      NumericVector seg, Eigen::MatrixXd W, Eigen::MatrixXd V,
                      NumericMatrix lambda1, NumericVector lambda2,
                      NumericMatrix alpha1, double alpha2,
                      bool multiply = false){
  int n = X.rows(), L = seg.size();

  NumericVector start_X(L), start_W(L), LL(L + 1);

  Eigen::VectorXd lambda21 = as<Eigen::Map<Eigen::VectorXd>>(lambda2);

  LL = LL_obs(X, Z, y, seg, W, V);

  NumericVector pen = Pen_L(seg, W, V, lambda1, lambda2, n,
                            alpha1, alpha2, multiply);
  return LL - pen;
}


//' Cpp function: Complete likelihood function of mixture distribution
//' @param X a data matrix of input in experts network.
//' @param Z a data matrix of input in gating network.
//' @param y a vector of response.
//' @param W parameters in experts network.
//' @param V parameters in gating network.
//' @return A matrix of complete likelihood by
//' compLikeli = pi * Likelihood = pi*(1/1+exp(-X*W))^y
//' pi = exp(Z * V_i)/sum(Z * V)
// [[Rcpp::export]]
Eigen::MatrixXd compLikeli(Eigen::MatrixXd X, Eigen::MatrixXd Z,
                            NumericVector y, Eigen::MatrixXd W,
                            Eigen::MatrixXd V){
  int n = X.rows(), p = X.cols(), K = V.cols(), q = Z.cols();
  Eigen::MatrixXd X1(n, p+1), Z1(n, q+1), ones = MatrixXd::Ones(n,1);
  X1 << ones, X; Z1 << ones, Z;
  Eigen::VectorXd prob(n);
  Eigen::MatrixXd latent_prob = calcProb(Z1, V);
  NumericMatrix sub_Likeli(n,K);
  NumericVector residual(n);
  residual.fill(1);
  for(int k=0; k<K; k++){
    sub_Likeli.column(k) = LBinom(X1, W.col(k), y, residual, false);
  }
  Eigen::MatrixXd sub_Likeli1 = as<Eigen::Map<Eigen::MatrixXd>>(sub_Likeli);
  Eigen::MatrixXd LL_nn;
  LL_nn = (latent_prob.array()*sub_Likeli1.array()).matrix();

  return LL_nn;
}

NumericVector LL_complete0(Eigen::MatrixXd X, Eigen::MatrixXd Z,
                           NumericVector y, Eigen::MatrixXd W,
                           Eigen::MatrixXd V){
  int n = X.rows(), K = V.cols();
  Eigen::MatrixXd LL_z, LL_n, LL_n1, LL_s, onekk = MatrixXd::Ones(K,K);
  Eigen::VectorXd Likeli(n), onek = MatrixXd::Ones(K,1);
  NumericVector prob1(n), Likelihood(n);
  LL_n = compLikeli(X, Z, y, W, V);
  LL_s = LL_n *onekk;
  LL_z = LL_n.array() / LL_s.array();
  LL_n1 = clippingMX(LL_n, 1e-9, 1);
  Likeli = (LL_n1.array().log() * LL_z.array()).matrix() * onek;
  Likelihood = wrap(Likeli);

  return Likelihood;
}

NumericVector LL_complete(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericMatrix y,
                          NumericVector seg, Eigen::MatrixXd W,
                          Eigen::MatrixXd V){
  int n = X.rows(), L = seg.size(), K = V.cols();

  NumericVector start_X(L), end_X(L), start_W(L), end_W(L), LL(L + 1);
  double LL_s = 0, LL_temp = 0;
  for(int l=0; l<L;l++){
    if(l == 0){
      start_X[l] = 0;
      start_W[l] = 0;
      end_X[l] = seg[l];
      end_W[l] = seg[l] + 1;
    }else{
      start_X[l] = end_X[l - 1];
      start_W[l] = end_W[l - 1];
      end_X[l]  = start_X[l] + seg[l];
      end_W[l]  = start_W[l] + seg[l] + 1;
    }
    LL_temp = sum(LL_complete0(X.block(0, start_X[l], n, seg[l]), Z,
                               y.column(l),  W.block(start_W[l], 0,
                               seg[l] + 1, K), V));

    LL[l] = LL_temp;
    LL_s += LL_temp;
  }
  LL[L] = LL_s;
  return LL;
}

NumericVector PLL_complete(Eigen::MatrixXd X, Eigen::MatrixXd Z,
                           NumericMatrix y, NumericVector seg,
                           Eigen::MatrixXd W, Eigen::MatrixXd V,
                           NumericMatrix lambda1, NumericVector lambda2,
                           NumericMatrix alpha1, double alpha2,
                           bool multiply = false){
  int n = X.rows(), L = seg.size();

  NumericVector start_X(L), start_W(L), LL(L + 1);

  Eigen::VectorXd lambda21 =  as<Eigen::Map<Eigen::VectorXd>>(lambda2);
  LL = LL_complete(X, Z, y, seg, W, V);

  NumericVector pen = Pen_L(seg, W, V, lambda1, lambda2,
                            n, alpha1, alpha2, multiply);
  return LL - pen;
}

//' Cpp function: calculate the four type of log likelihood function
//' (observed likelihood, penalized likelihood,
//' complete likelihood, penalized complete likelihood)
//' in all levels.
//' @param X an aggregated data matrix (n*P) of input in experts network
//' (rbind of all levels input).
//' @param Z a data matrix (n*q) of input in gating network.
//' @param y a vector of response (n).
//' @param seg an integer vector of length in each level (L).
//' @param W aggregated parameters in experts network (cbind of all levels)
//' ((P+L)*K).
//' @param V parameters in gating network((q+1)*K).
//' @param lambda1 aggregated penalty lambda parameters in experts network(P*K).
//' @param lambda2 penalty parameters lambda in gating network(q*1).
//' @param alpha1 penalty parameters alpha in experts network(L*K).
//' @param alpha2 penalty parameter alpha in gating network(1).
//' @param multiply a bool variable indicate whether lambda parameters
//'  times sample size.
//' @return A matrix of four types of likelihood function.
//'
// [[Rcpp::export]]
NumericMatrix calcLL(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericMatrix y,
                     NumericVector seg, Eigen::MatrixXd W, Eigen::MatrixXd V,
                     NumericMatrix lambda1, NumericVector lambda2,
                     NumericMatrix alpha1, double alpha2,
                     bool multiply = false){
  int L = seg.size();
  NumericMatrix L_a((L+1), 4);
  L_a(_,0) = LL_obs(X, Z, y, seg, W, V);
  L_a(_,1) = PLL_obs(X, Z, y, seg, W, V,
            lambda1, lambda2, alpha1, alpha2, multiply);
  L_a(_,2) = LL_complete(X, Z, y, seg, W, V);
  L_a(_,3) = PLL_complete(X, Z, y, seg, W, V,
            lambda1, lambda2, alpha1, alpha2, multiply);

  return L_a;
}
/* -------------------------------------------
 ---------- Elnet Function ------------
 ------------------------------------------- */

Eigen::VectorXd Elnet_step(Eigen::MatrixXd X, Eigen::VectorXd y,
                           Eigen::VectorXd R, Eigen::VectorXd beta_old,
                           Eigen::VectorXd lambda, double alpha){


  int p = beta_old.size();
  double beta_temp, Xi_norm, R_s = R.sum();
  Eigen::VectorXd Xi, beta_i = beta_old, beta_new = beta_old;
  Eigen::MatrixXd R_mat = R.asDiagonal();


  for(int i =0; i<p; i++){
    Xi = X.col(i);
    beta_i(i) = 0;
    Xi_norm  = (R_mat * Xi).dot(Xi);
    if(i == 0){
      beta_temp = (R_mat*(y - X * beta_i)).sum() / R_s ;
    }else{
      beta_temp = ((y - X * beta_i).dot(R_mat*Xi));
      beta_temp = softThresh(beta_temp, lambda(i - 1)*alpha);
      beta_temp = double(beta_temp) / (Xi_norm + lambda(i - 1)*(1 - alpha));
    }
    beta_i(i) = beta_temp;
  }
  return beta_i;
}

Eigen::VectorXd Elnet1(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd R,
                       Eigen::VectorXd lambda, double alpha,
                       int itmax = 1e3, double eps_var = 1e-5){

  int n = X.rows(), p = X.cols();
  Rcpp::NumericVector y_1 = wrap(y);
  Eigen::VectorXd Beta0;
  Eigen::MatrixXd X1(n, p +1);
  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n,1);
  X1 << ones, X;


  Beta0=Eigen::VectorXd::Zero(p + 1);
  Beta0(0) = mean(y_1);


  Eigen::VectorXd Beta = Beta0;
  lambda = lambda * n;

  if(var(y_1) < eps_var){
    return Beta0;
  }else{

    int it = 0;
    bool flag = true;
    while(flag){

      Beta = Elnet_step(X1, y, R, Beta0, lambda, alpha);

      bool cond1 = (double(1)/n * (Beta - Beta0).dot(Beta - Beta0) < 1e-9);
      bool cond2 = (it > itmax);

      if(cond1 || cond2){
        flag = false;
      }

      Beta0 = Beta;
      it = it + 1;
    }
    return Beta;
  }
}

Eigen::VectorXd Elnet0(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd R,
                       Eigen::VectorXd lambda, double alpha){
  //Rcpp::List myglmnet(Eigen::MatrixXd x, Eigen::VectorXd y, double lambda){

  int n = X.rows();
  int p = X.cols();

  double lambda1;
  lambda1 = lambda.sum() * double(n) / (double(p) * R.sum());

  Eigen::VectorXd pen = lambda * p/ lambda.sum();

  Rcpp::Environment pkg_glmnet = Environment::namespace_env("glmnet");
  Rcpp::Function glmnet = pkg_glmnet["glmnet"];
  Rcpp::List glmres = glmnet(Rcpp::Named("x",X), Rcpp::Named("y",y),
                             Rcpp::Named("family", "gaussian"),
                             Rcpp::Named("lambda", lambda1),
                             Rcpp::Named("alpha", alpha),
                             Rcpp::Named("standardize", false),
                             Rcpp::Named("penalty.factor",pen),
                             Rcpp::Named("intercept", true),
                             Rcpp::Named("weights", R));

  Eigen::SparseMatrix<double> glmcoef0 = glmres["beta"];
  Eigen::VectorXd a0 = glmres["a0"];
  Eigen::VectorXd beta(p + 1);
  beta.head(1) = a0;
  beta.tail(p) = glmcoef0;

  return beta;
}

Eigen::VectorXd Elnet(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd R,
                      Eigen::VectorXd lambda, double alpha){

  Eigen::VectorXd beta;

  try{
    beta = Elnet0(X, y, R, lambda, alpha);
  }catch(...){
    beta = Elnet1(X, y, R, lambda, alpha);
  }
  return beta;

}

double Elobj(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd R,
             Eigen::VectorXd beta, Eigen::VectorXd lambda,
             double alpha, bool multiply){

  int n = X.rows(), p = X.cols();
  Eigen::MatrixXd R_mat = R.asDiagonal();
  Eigen::VectorXd B(p);
  double Loss = 0, pen = 0;
  Eigen::MatrixXd X1(n, p + 1);
  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n,1);
  X1 << ones, X;
  B = beta.tail(p);

  if(multiply){
    lambda = lambda * n;
  }

  Loss = 0.5*(y - X1 * beta).dot(R_mat * (y - X1 * beta));

  for(int i =0; i < p; i++){
    pen += lambda(i) * (0.5*(1 - alpha)*pow(B(i),2) + alpha * abs(B(i)));
  }
  return Loss + pen;
}


/* -------------------------------------------
 ---------- sMulti Function ------------
 ------------------------------------------- */

double obj_logit(Eigen::MatrixXd X, Eigen::VectorXd beta, NumericVector y,
                 NumericVector R, Eigen::VectorXd lambda, double alpha,
                 bool multiply, NumericVector residual_k){

  int n = X.rows(), p = X.cols();
  Eigen::MatrixXd X1(n, p + 1);
  Eigen::VectorXd B(p);
  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n,1);
  double Loss = 0, pen = 0;
  X1 << ones, X;

  if(multiply){
    lambda = lambda * n;
  }
  B = beta.tail(p);


  Loss = sum(R * LBinom(X1, beta, y, residual_k, true));

  for(int i =0; i<p; i++){
    pen += lambda(i) * (0.5*(1 - alpha)*pow(B(i),2) + alpha * abs(B(i)));
  }
  return -Loss + pen;

}

Eigen::VectorXd BtrLogit(Eigen::MatrixXd X, NumericVector y, NumericVector R,
                         Eigen::VectorXd lambda, double alpha,
                         Eigen::VectorXd beta_old, Eigen::VectorXd beta_new,
                         NumericVector residual_k, double shrink = 0.8){


  Eigen::VectorXd d_x = beta_new - beta_old;

  double f_old,f_new;
  f_old = obj_logit(X, beta_old, y, R, lambda, alpha, true, residual_k);
  f_new = obj_logit(X, beta_new, y, R, lambda, alpha, true, residual_k);

  NumericVector f_val = NumericVector::create(f_old, f_new);

  LogicalVector inf_val, na_val, na_new;
  bool cond1, cond2;
  inf_val = is_infinite(f_val);
  na_val = is_nan(f_val);

  bool cond_old = is_true(any(inf_val)) || is_true(any(na_val));
  int it = 0;
  double t = 1;
  bool flag = true;

  if(cond_old != 0){
    return beta_old;
  }

  while(flag){
    beta_new = beta_old + t * d_x;

    f_new = obj_logit(X, beta_new, y, R, lambda, alpha, true, residual_k);
    f_val(1) = f_new;

    na_new = is_nan(f_val);
    if(is_true(any(na_new))){
      return beta_old;
    }

    cond1 = ((f_new - f_old) < (1e-9));

    cond2 = (it > 100);

    if(cond1 || cond2){

      flag = false;
    }
    t = shrink * t;
    it = it + 1;
  }

  return beta_new;

}
Eigen::VectorXd sMulti_single(Eigen::MatrixXd X, Eigen::VectorXd y,
                              Eigen::VectorXd beta0, NumericVector R,
                              Eigen::VectorXd lambda, NumericVector residual_k,
                              double alpha, double beta_max,
                              int itmax = 1e3, bool btr = true){

  int n = X.rows(), p = X.cols();
  Eigen::VectorXd beta = beta0;
  Eigen::MatrixXd X1(n, p+1);
  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n,1);
  Eigen::VectorXd z0, weight0;
  X1 << ones, X;
  int it = 0;
  bool cond1, cond2, flag = true;
  Eigen::VectorXd beta_btr, odd(n);
  NumericVector odd1, pi_tid, w_tid, weight, z, beta1, y1 = wrap(y);

  while(flag){

    odd = X1 * beta0;
    odd1 = wrap(odd);
    pi_tid = 1 / (1 + residual_k * Rcpp::exp(-odd1));
    pi_tid = checking(pi_tid);

    w_tid = pi_tid * (1 - pi_tid);
    w_tid = clippingV(w_tid, 1e-9, 1 - 1e-9);

    weight = w_tid * R;

    z = odd1 + (y1 - pi_tid) / w_tid;
    z0 = as<Eigen::Map<Eigen::VectorXd>>(z);
    weight0 = as<Eigen::Map<Eigen::VectorXd>>(weight);
    beta = Elnet(X, z0, weight0, lambda, alpha);
    beta_btr = clippingVX(beta, -beta_max, beta_max);
    if(btr){

      beta_btr = BtrLogit(X, y1, R, lambda, alpha, beta0,
                          beta_btr, residual_k, 0.8);
    }

    cond1 = ((beta_btr - beta0).array().abs().sum() < 1e-4);
    cond2 = (it > itmax);
    if(cond1 || cond2){
      flag = false;
    }
    beta0 = beta_btr;
    it = it + 1;
  }
  return beta0;
}

Eigen::VectorXd sLogit1(Eigen::MatrixXd X, Eigen::VectorXd y,
                        Eigen::VectorXd lambda, NumericVector R,
                        double alpha, double beta_max,
                        int itmax = 1e3, bool btr = true){

  int p = X.cols(), n = X.rows();
  double pi_temp;
  Eigen::VectorXd beta0(p + 1), beta(p + 1);
  NumericVector residual_k(n);
  residual_k.fill(1);
  pi_temp = y.array().mean();
  pi_temp = clipping(pi_temp, 1e-5, 1 - 1e-5);
  beta0 = Eigen::MatrixXd::Zero(p + 1, 1);
  beta0(0) = log(pi_temp/(1 - pi_temp));
  if(abs(beta0(0)) >= -log(1e-5/1-(1e-5))){
    return beta0;
  }else{
    beta = sMulti_single(X, y, beta0, R, lambda, residual_k,
                         alpha, beta_max,  itmax, btr);
    return beta;
  }
}

Eigen::MatrixXd sMulti_step(Eigen::MatrixXd X, Eigen::MatrixXd y,
                            Eigen::MatrixXd Beta0,
                            Eigen::VectorXd lambda, NumericVector R,
                            double alpha, double beta_max,
                            int itmax = 1e3, bool btr = true){

  int n = X.rows(), p = X.cols(), p1 = Beta0.rows(), K = y.cols();

  Eigen::MatrixXd Beta_k(p1, K - 1), Beta(p1, K);
  Eigen::MatrixXd X1(n, p +1);
  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n,1),
    oneK1 = Eigen::MatrixXd::Ones(K - 1,1);
  Eigen::VectorXd residual_k, Betak, yk, beta;
  NumericVector residualk;
  X1 << ones, X;
  Beta = Beta0;

  for(int k = 0; k < K; k++){

    if(abs(Beta0(0, k)) <= -log((1e-5)/(1-(1e-5)))){
      Beta_k = Exccol(Beta, k);
      residual_k = (X1 * Beta_k).array().exp().matrix()*oneK1;
      residualk = wrap(residual_k);
      Betak = Beta.col(k);
      yk = y.col(k);
      beta = sMulti_single(X, yk, Betak, R, lambda,
                           residualk, alpha, beta_max,
                           itmax, btr);
      beta = clippingVX(beta, -beta_max, beta_max);
      Beta.col(k) = beta;
    }

  }
  return Beta;

}

Eigen::MatrixXd sMulti1(Eigen::MatrixXd X, Eigen::MatrixXd y,
                        Eigen::VectorXd lambda, NumericVector R,
                        double alpha, double beta_max, int itmax1 = 1e2,
                        int itmax2 = 1e3, bool btr = true){

  int n = X.rows(), p = X.cols(), K0=y.cols();
  Eigen::MatrixXd X1(n, p +1);
  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n,1),
    onep1 = Eigen::MatrixXd::Ones(p + 1,1);
  Eigen::MatrixXd Beta0, Beta;
  Eigen::MatrixXd y1(n, 2);
  NumericVector beta1;

  double pi_temp;
  Eigen::ArrayXd pi_temp1;
  int K = K0;
  if(K0 == 1){
    y1 << ones - y, y;
    K = 2;
  }
  LogicalVector early_stop(K - 1);

  X1 << ones, X;
  Beta0 = Eigen::MatrixXd::Zero(p + 1, K);
  Beta = Eigen::MatrixXd::Zero(p + 1, K);

  for(int k =0; k < K; k ++){

    if(K0 == 1){
      pi_temp = y1.col(k).sum()/n;
    }else{
      pi_temp = y.col(k).sum()/n;
    }
    early_stop[k] = (pi_temp < 1e-5) || (pi_temp > 1-1e-5);
    pi_temp = clipping(pi_temp, 1e-5, 1 - 1e-5);
    Beta0(0,k) =log(pi_temp/(1 - pi_temp));
  }

  if(is_true(any(early_stop))){
    return Beta0;
  }


  int it = 0;
  bool cond1 = false, cond2 = false, flag = true;
  while (flag){

    if(K0 == 1){
      Beta = sMulti_step(X, y1, Beta0, lambda, R, alpha, beta_max, itmax2, btr);
    }else{
      Beta = sMulti_step(X, y, Beta0, lambda, R, alpha, beta_max, itmax2, btr);
    }

    for(int k =0; k < p + 1; k++){
      beta1 = wrap(Beta.row(k));
      if(k == 0){
        beta1 = beta1 - mean(beta1);
      }else{
        beta1 = beta1 - median(beta1);
      }
      Beta.row(k) = as<Eigen::Map<Eigen::MatrixXd>>(beta1);
    }

    cond1 = ((Beta - Beta0).array().abs().sum() < 1e-9);
    cond2 = (it > itmax1);
    if( cond1 || cond2){
        flag = false;
    }
    it = it + 1;
    Beta0 = Beta;
  }

  return Beta;
}

Eigen::VectorXd sLogit0(Eigen::MatrixXd X, Eigen::VectorXd y,
                        Eigen::VectorXd lambda, NumericVector R,
                        double alpha, double beta_max){

  int n = X.rows(), p = X.cols();

  double lambda1;
  Eigen::MatrixXd zero = Eigen::MatrixXd::Zero(p + 1, 1);
  lambda1 = lambda.sum() * double(n) / (double(p) * sum(R));

  Eigen::VectorXd pen = lambda * p/ lambda.sum();

  Rcpp::Environment pkg_glmnet = Environment::namespace_env("glmnet");
  Rcpp::Function glmnet = pkg_glmnet["glmnet"];
  Rcpp::List glmres = glmnet(Rcpp::Named("x",X), Rcpp::Named("y",y),
                             Rcpp::Named("family", "binomial"),
                             Rcpp::Named("lambda", lambda1),
                             Rcpp::Named("alpha", alpha),
                             Rcpp::Named("standardize", false),
                             Rcpp::Named("penalty.factor",pen),
                             Rcpp::Named("intercept", true),
                             Rcpp::Named("weights", R));

  Eigen::SparseMatrix<double> glmcoef0 = glmres["beta"];
  Eigen::VectorXd a0 = glmres["a0"];
  Eigen::VectorXd beta(p + 1);
  beta.head(1) = a0;
  beta.tail(p) = glmcoef0;
  beta = clippingVX(beta, -beta_max, beta_max);

  return beta;
}

Eigen::VectorXd sLogit(Eigen::MatrixXd X, Eigen::VectorXd y,
                       Eigen::VectorXd lambda,
                       NumericVector R, double alpha, double beta_max){
  Eigen::VectorXd beta;

  try{
    beta = sLogit0(X, y, lambda, R, alpha, beta_max);
  }catch(...){
    beta = sLogit1(X, y, lambda, R, alpha, beta_max);
  }
  return beta;
}

Eigen::MatrixXd sMulti0(Eigen::MatrixXd X, Eigen::MatrixXd y, Eigen::VectorXd lambda,
                        NumericVector R, double alpha, double beta_max){

  int n = X.rows(), p = X.cols(), K0 = y.cols();
  Eigen::MatrixXd y1(n, 2);
  Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(n,1);
  int K = K0;

  if(K0 == 1){
    y1 << ones - y, y;
    K = 2;
  }
  Eigen::MatrixXd Beta(p +1, K);

  double lambda1;
  lambda1 = lambda.sum() * double(n) / (double(p) * sum(R));

  Eigen::VectorXd pen = lambda * p/ lambda.sum();

  Rcpp::Environment pkg_glmnet = Environment::namespace_env("glmnet");
  Rcpp::Function glmnet = pkg_glmnet["glmnet"];
  Rcpp::List glmres = glmnet(Rcpp::Named("x",X), Rcpp::Named("y",y),
                             Rcpp::Named("family", "multinomial"),
                             Rcpp::Named("lambda", lambda1),
                             Rcpp::Named("alpha", alpha),
                             Rcpp::Named("standardize", false),
                             Rcpp::Named("penalty.factor",pen),
                             Rcpp::Named("intercept", true),
                             Rcpp::Named("weights", R));

  NumericVector beta(p+1), beta1(p);
  NumericVector a0 = glmres["a0"];
  List BETA = glmres["beta"];
  Eigen::SparseMatrix<double> glmcoef0;
  Eigen::MatrixXd glmcoef;

  for(int k =0; k<K;k++){
    beta[0] = a0[0];
    glmcoef0 = BETA[k];
    glmcoef = glmcoef0;
    beta1 = wrap(glmcoef);
    beta[seq(1,p)] = beta1;
    Beta.col(k) = as<Eigen::Map<Eigen::MatrixXd>>(beta);
  }
  return Beta;
}

//' Cpp function: fitting sparse multinomial regression.
//' @param X input matrix for sparse multinomial regression(n*p).
//' @param y response data for sparse multinomial regression(n*K).
//' @param lambda A vector of penalty parameters lambda (p*1).
//' @param alpha A number of penalty parameters alpha.
//' @param R A vector of weighted parameters for each sample(n*1).
//' @param beta_max A number of maximal number of fitted coefficients.
//' @return A matrix of fitted coefficients.
// [[Rcpp::export]]
Eigen::MatrixXd sMulti(Eigen::MatrixXd X, Eigen::MatrixXd y,
                       Eigen::VectorXd lambda, NumericVector R,
                       double alpha, double beta_max){
  Eigen::MatrixXd Beta;

  try{
    Beta = sMulti0(X, y, lambda, R, alpha, beta_max);
  }catch(...){
    Beta = sMulti1(X, y, lambda, R, alpha, beta_max);
  }
  return Beta;
}

/* -------------------------------------------
 ---------- EM Function ------------
 ------------------------------------------- */

Eigen::MatrixXd NEMoE_CEM_W(Eigen::MatrixXd X, Eigen::MatrixXd Z,
                            NumericVector y, NumericMatrix r_i,
                            NumericMatrix lambda1, NumericVector alpha1,
                            double beta_max, int itmax = 1e3,
                            double gamma = 1){

  int K = r_i.cols(), p = X.cols(), q = Z.cols();
  IntegerMatrix r_hard = hard_prob(r_i);
  Eigen::MatrixXd X_temp, W = Eigen::MatrixXd::Zero(p + 1,K),
    V= Eigen::MatrixXd::Zero(q + 1, K);
  Eigen::MatrixXd lambda10 = as<Eigen::Map<Eigen::MatrixXd>>(lambda1);
  NumericVector y_temp;
  Eigen::VectorXd y_temp0, W_k;

  for(int k=0;k<K;k++){

    IntegerVector Pos = pos1(r_hard.column(k));

    X_temp = RowSelX(X, Pos);
    y_temp = y[Pos];
    y_temp0 = as<Eigen::Map<Eigen::VectorXd>>(y_temp);

    NumericVector R(Pos.size());
    R.fill(1);
    W_k = sLogit(X_temp, y_temp0, lambda10.col(k), R, alpha1[k], beta_max);

    W.col(k) = W_k;
  }

  return W;
}

Eigen::MatrixXd NEMoE_BEM_W(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericVector y,
                            NumericMatrix r_i, NumericMatrix lambda1,
                            NumericVector alpha1, double beta_max,
                            int itmax = 5, double gamma = 1){
  int K = r_i.cols(), p = X.cols(), q = Z.cols();
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(p + 1,K),
    V= Eigen::MatrixXd::Zero(q + 1, K), W_k = Eigen::MatrixXd::Zero(p + 1,1);
  Eigen::MatrixXd lambda10 = as<Eigen::Map<Eigen::MatrixXd>>(lambda1);
  Eigen::VectorXd y0 = as<Eigen::Map<Eigen::VectorXd>>(y);

  for(int k=0;k<K;k++){

    W_k = sLogit(X, y0, lambda10.col(k), r_i.column(k), alpha1[k], beta_max);

    W.col(k) = W_k;
  }

  return W;
}

Eigen::MatrixXd NEMoE_SEM_W(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericVector y,
                            NumericMatrix r_i, NumericMatrix lambda1,
                            NumericVector alpha1, double beta_max,
                            int itmax = 5, double gamma = 1){

  NumericMatrix r_sample = rSample(r_i);
  Eigen::MatrixXd W = NEMoE_CEM_W(X, Z, y, r_sample, lambda1, alpha1,
                                  beta_max, itmax,  gamma);
  return W;

}

Eigen::MatrixXd NEMoE_GEM_W(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericVector y,
                            NumericMatrix r_i, NumericMatrix lambda1,
                            NumericVector alpha1, double beta_max,
                            int itmax = 5, double gamma = 1){

  int K = r_i.cols(), p = X.cols(), q = Z.cols();
  Eigen::MatrixXd W(p + 1,K), V(q + 1, K);
  Eigen::MatrixXd lambda10 = as<Eigen::Map<Eigen::MatrixXd>>(lambda1);
  Eigen::VectorXd W_k, y0 = as<Eigen::Map<Eigen::VectorXd>>(y);

  // #pragma omp parallel for
  for(int k=0;k<K;k++){

    W_k = sLogit1(X, y0, lambda10.col(k), r_i.column(k), alpha1[k],
                                  beta_max, itmax);

    W.col(k) = W_k;
  }
  return W;
}

Eigen::MatrixXd NEMoE_SAEM_W(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericVector y,
                             NumericMatrix r_i, NumericMatrix lambda1,
                             NumericVector alpha1, double beta_max,
                             int itmax = 5, double gamma = 1){

  NumericMatrix r_i_copy = clone(r_i);
  NumericMatrix r_1 = rSample(r_i_copy);

  Eigen::MatrixXd r_10 = as<Eigen::Map<Eigen::MatrixXd>>(r_1);
  Eigen::MatrixXd r_i0 = as<Eigen::Map<Eigen::MatrixXd>>(r_1);
  Eigen::MatrixXd r_saem0 = r_10 * (1 - gamma) + r_i0 * gamma;
  NumericMatrix r_saem = wrap(r_saem0);

  Eigen::MatrixXd W = NEMoE_BEM_W(X, Z, y, r_i, lambda1,
                                  alpha1, beta_max, itmax, gamma);

  return W;
}

Eigen::MatrixXd NEMoE_EM1(Eigen::MatrixXd X, NumericVector seg,
                          Eigen::MatrixXd Z, NumericMatrix y, NumericMatrix r_i,
                          NumericMatrix lambda1, NumericVector lambda2,
                          NumericMatrix alpha1, double alpha2,
                          double beta_max, int EM_opt, double gamma){

  int n = X.rows(), P = X.cols(), L = seg.size(), q = Z.cols(), K = r_i.cols();
  Eigen::MatrixXd W(P+L,K), V(q+1,K), r_i_s0;
  Eigen::VectorXd V_temp(q+1);
  NumericVector R0(n);
  R0.fill(1);
  Eigen::MatrixXd r_i1, lambda20 = as<Eigen::Map<Eigen::MatrixXd>>(lambda2),
    MoE_res(P + L + q + 1, K);
  NumericVector start_X(L), end_X(L), start_W(L), end_W(L);
  NumericMatrix r_i_l, lambda1_temp, r_i_s;
  for(int l=0; l<L;l++){
    if(l == 0){
      start_X[l] = 0;
      start_W[l] = 0;
      end_X[l] = seg[l];
      end_W[l] = seg[l] + 1;
    }else{
      start_X[l] = end_X[l - 1];
      start_W[l] = end_W[l - 1];
      end_X[l]  = start_X[l] + seg[l];
      end_W[l]  = start_W[l] + seg[l] + 1;
    }
  }

  for(int l=0; l<L; l++){
    r_i_l = RowSel(r_i, seq(l*n, (l + 1)*n - 1));
    lambda1_temp = lambda1(seq(start_X[l], end_X[l] - 1),_);
    if(EM_opt == 1){
      W.block(start_W[l], 0, seg[l] + 1, K) =
        NEMoE_CEM_W(X.block(0, start_X[l], n, seg[l]),
                            Z, y.column(l), r_i_l, lambda1_temp,
                            alpha1.row(l), beta_max);
    }else if(EM_opt == 2){
      W.block(start_W[l], 0, seg[l] + 1, K) =
        NEMoE_SEM_W(X.block(0, start_X[l], n, seg[l]),
                    Z, y.column(l), r_i_l, lambda1_temp,
                    alpha1.row(l), beta_max);
    }else if(EM_opt == 3){
      W.block(start_W[l], 0, seg[l] + 1, K) =
        NEMoE_SAEM_W(X.block(0, start_X[l], n, seg[l]),
                     Z, y.column(l), r_i_l, lambda1_temp,
                     alpha1.row(l), beta_max, 5, gamma);
    }else if(EM_opt == 4){
      W.block(start_W[l], 0, seg[l] + 1, K) =
        NEMoE_GEM_W(X.block(0, start_X[l], n, seg[l]),
                    Z, y.column(l), r_i_l, lambda1_temp,
                    alpha1.row(l), beta_max);
    }else{
      W.block(start_W[l], 0, seg[l] + 1, K) =
        NEMoE_BEM_W(X.block(0, start_X[l], n, seg[l]),
                    Z, y.column(l), r_i_l, lambda1_temp,
                    alpha1.row(l), beta_max);
    }
  }
  r_i_s = R_sum(r_i, L);
  r_i_s0 = as<Eigen::Map<Eigen::MatrixXd>>(r_i_s);

  if(K == 2){
    V_temp = sLogit1(Z, r_i_s0.col(1), lambda20, R0, alpha2, beta_max);
    V.col(0) = (-V_temp/2);
    V.col(1) = V_temp/2;
  }else{
    V = sMulti1(Z, r_i_s0, lambda20, R0, alpha2, beta_max);
  }

  MoE_res << W,V;

  return MoE_res;
}

double obj_W(Eigen::MatrixXd X, NumericVector y,
             NumericVector r_i, NumericVector lambda1,
             double alpha1, Eigen::VectorXd W,
             Eigen::VectorXd pen_fac){

  int n = X.rows(), p = X.cols();
  Eigen::MatrixXd X1(n, p+1), ones = MatrixXd::Ones(n,1);
  Eigen::VectorXd lambda11;
  lambda11 = as<Eigen::Map<Eigen::VectorXd>>(lambda1);
  X1 << ones, X;

  Eigen::VectorXd lambda11_fac = (lambda11.array() * pen_fac.array()).matrix();
  NumericVector residual_k(n);
  residual_k.fill(1);
  double PLL = -obj_logit(X, W, y, r_i, lambda11_fac, alpha1, true, residual_k);

  return PLL;

}

Eigen::MatrixXd calcPen_fac(Eigen::MatrixXd X, Eigen::MatrixXd Z,
                            NumericVector y, Eigen::MatrixXd W,
                            Eigen::MatrixXd V, bool equal_var = false){

  int n = X.rows(), p = X.cols(), K = V.cols(), q= Z.cols();
  Eigen::MatrixXd mu(p,K), sgm(p,K), pen_fac(p,K), r_temp(p,K), r_i,
  onepl = MatrixXd::Ones(p,1), ones = MatrixXd::Ones(n,1);
  Eigen::MatrixXd Z1(n,q+1), pi(n,K);
  Z1 << ones,Z;

  if(equal_var){

    pi = calcProb(Z1,V);
    pen_fac = onepl * ones.transpose() * pi/n;

  }else{
//    r_i = compLikeli(X, Z, y, W, V);
    r_i = calcProb(Z1,V);
    r_temp = onepl * ones.transpose() * r_i;
    mu = (((X.transpose()* r_i).array())/(r_temp.array())).matrix();
    sgm = (((X.array().square().matrix().transpose() * r_i).array() -
      ((mu.array().square()) * (r_temp.array())).sqrt())
             /(r_temp.array())).matrix();
    pen_fac = ((((onepl * ones.transpose() * r_i).array())/n)
                 *sgm.array()).matrix();

  }

  return(pen_fac);
}

Eigen::MatrixXd calcPen_fac1(Eigen::MatrixXd X, NumericVector seg,
                             Eigen::MatrixXd Z, NumericMatrix y,
                             Eigen::MatrixXd W, Eigen::MatrixXd V,
                             bool adapt){

  int L = seg.size(), P = X.cols(), K = V.cols(), n = X.rows();
  Eigen::MatrixXd pen_fac(P,K), onepK= MatrixXd::Ones(P, K);
  NumericVector start_X(L), end_X(L), start_W(L), end_W(L);

  for(int l=0; l<L;l++){
    if(l == 0){
      start_X[l] = 0;
      start_W[l] = 0;
      end_X[l] = seg[l];
      end_W[l] = seg[l] + 1;
    }else{
      start_X[l] = end_X[l - 1];
      start_W[l] = end_W[l - 1];
      end_X[l]  = start_X[l] + seg[l];
      end_W[l]  = start_W[l] + seg[l] + 1;
    }
  }
  if(adapt){
    for(int l = 0; l<L; l++){
      pen_fac.block(start_X[l],0, seg[l],K) = calcPen_fac(
        X.block(0, start_X[l], n, seg[l]), Z, y.column(l),
        W.block(start_W[l], 0, seg[l] + 1, K), V);
    }

  }else{
    pen_fac = double(1)/double(K) * onepK;
  }

  return(pen_fac);
}

Eigen::MatrixXd NEMoE_btrW(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericVector y,
                           NumericMatrix r_i, NumericMatrix lambda1,
                           NumericVector alpha1, Eigen::MatrixXd W,
                           Eigen::MatrixXd W0, Eigen::MatrixXd V,
                           Eigen::MatrixXd pen_fac, Eigen::MatrixXd pen_fac0,
                           bool adapt, double shrink = 0.8){

  int K = W.cols();
  double f_old, f_new, t;
  NumericVector f_val;
  LogicalVector inf_val, na_val, na_new;
  bool cond1, cond2, cond_old, flag;
  int it;
  Eigen::VectorXd W0_temp, W_temp, d_W_temp;
  Eigen::MatrixXd pen_fac1;
  for(int k=0; k<K; k++){

    W0_temp = W0.col(k);
    W_temp = W.col(k);
    d_W_temp = W_temp - W0_temp;

    f_old = obj_W(X, y, r_i.column(k), lambda1.column(k),
                  alpha1[k], W0_temp, pen_fac0.col(k));
    f_new = obj_W(X, y, r_i.column(k), lambda1.column(k),
                  alpha1[k], W_temp, pen_fac0.col(k));
    f_val = {f_old, f_new};

    inf_val = is_infinite(f_val);
    na_val = is_nan(f_val);
    cond_old = is_true(any(inf_val)) || is_true(any(na_val));

    it = 0;
    t = 1;

    W_temp = W0.col(k);
    if(cond_old != 0){
      W.col(k) = W0_temp;
    }else{
      flag = true;
      while(flag){
        W_temp = W0_temp + t*d_W_temp;

//        if(adapt){
//          pen_fac1 = calcPen_fac(X, Z, y, W, V);
//        }else{
//          pen_fac1 = pen_fac;
//        }
        f_new = obj_W(X, y, r_i.column(k), lambda1.column(k),
                      alpha1[k], W_temp, pen_fac0.col(k));
        f_val(1) = f_new;
        na_new = is_nan(f_val);
        if(is_true(any(na_new))){
          W.col(k) = W0_temp;
          flag = false;
        }else{

          cond1 = ((f_new - f_old) > -(1e-9));

          cond2 = (it > 20);

          if(cond1 || cond2){

            flag = false;
          }else{
            W = W0;
          }
          t = shrink * t;
          it = it + 1;
        }
      }
    }
  }
  return W;
}

double obj_V(Eigen::MatrixXd Z, NumericMatrix r_i_s, NumericVector lambda2,
             double alpha2, Eigen::MatrixXd V){

  int n = Z.rows(), q = Z.cols();
  Eigen::MatrixXd Z1(n, q+1), ones = MatrixXd::Ones(n,1);
  Z1 << ones, Z;
  Eigen::MatrixXd r_i_s1, latent_prob = calcProb(Z1, V);
  r_i_s1 = as<Eigen::Map<Eigen::MatrixXd>>(r_i_s);
  double LL = (latent_prob.array() * log(r_i_s1.array())).sum();
  double penalty = pen_M(V, lambda2 * n, alpha2);
  return LL - penalty;
}

double obj_V1(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericMatrix y,
              NumericVector seg, NumericMatrix lambda1, NumericMatrix alpha1,
              NumericVector lambda2, double alpha2, Eigen::MatrixXd W,
              Eigen::MatrixXd V){

  int n = Z.rows(), L = seg.size(), P = X.cols(), K = V.cols();
  Eigen::MatrixXd ones = MatrixXd::Ones(n,1), pen_fac(P,K), lambda1X(P,K);
  NumericVector PLL;
  NumericMatrix lambda1_eff;
  double PLL0;
  pen_fac = calcPen_fac1(X, seg, Z, y, W, V, true);
  lambda1X = as<Eigen::Map<Eigen::MatrixXd>>(lambda1);
  lambda1_eff = wrap(pen_fac.array() * lambda1X.array());

  PLL = PLL_obs(X, Z, y, seg, W, V, lambda1_eff, lambda2, alpha1, alpha2, true);
  PLL0 = PLL[L];

  return PLL0;
}

Eigen::MatrixXd NEMoE_btrV(Eigen::MatrixXd X, Eigen::MatrixXd Z, NumericMatrix y,
                           NumericVector seg, NumericMatrix r_i,
                           NumericMatrix lambda1, NumericMatrix alpha1,
                           NumericVector lambda2, double alpha2, Eigen::MatrixXd W,
                           Eigen::MatrixXd V, Eigen::MatrixXd V0,
                           bool adapt, double PLL0, double shrink = 0.8){

  int L = seg.size();
  V = projMultiConstraint(V, alpha2);
  V0 = projMultiConstraint(V0, alpha2);
  Eigen::MatrixXd d_V = V - V0;
  double f_old, f_new;
  NumericMatrix r_i_s;
  r_i_s = R_sum(r_i, L);
  if(adapt){
    f_old = PLL0;
    f_new = obj_V1(X, Z, y, seg, lambda1, alpha1, lambda2, alpha2, W, V);
  }else{
    f_old = obj_V(Z, r_i_s, lambda2/double(L), alpha2, V0);
    f_new = obj_V(Z, r_i_s, lambda2/double(L), alpha2, V);
  }

  NumericVector f_val = NumericVector::create(f_old, f_new);

  LogicalVector inf_val, na_val, na_new;
  bool cond1, cond2;
  inf_val = is_infinite(f_val);
  na_val = is_nan(f_val);

  bool cond_old = is_true(any(inf_val)) || is_true(any(na_val));
  int it = 0;
  double t = 1;
  bool flag = true;

  if(cond_old != 0){
    return V0;
  }

  while(flag){
    V = V0 + t * d_V;
    V = projMultiConstraint(V, alpha2);
    if(adapt){
      f_new = obj_V1(X, Z, y, seg, lambda1, alpha1, lambda2, alpha2, W, V);
    }else{
      f_new = obj_V(Z, r_i_s, lambda2/double(L), alpha2, V);
    }
    f_val(1) = f_new;

    na_new = is_nan(f_val);
    if(is_true(any(na_new))){
      return V0;
    }

    cond1 = ((f_new - f_old) > -(1e-9));

    cond2 = (it > 20);

    if(cond1 || cond2){

      flag = false;
    }else{
      V = V0;
    }
    t = shrink * t;
    it = it + 1;
  }

  return V;
}

Eigen::MatrixXd NEMoE_step(Eigen::MatrixXd X, NumericVector seg,
                           Eigen::MatrixXd Z, NumericMatrix y,
                           NumericMatrix lambda1, NumericVector lambda2,
                           NumericMatrix alpha1, double alpha2,
                           Eigen::MatrixXd V_old, Eigen::MatrixXd W_old,
                           double beta_max, int EM_opt, Eigen::MatrixXd pen_fac0,
                           bool adapt, double gamma, bool btr, double PLL0){

  int n = X.rows(), P = X.cols(), L = seg.size(), q = Z.cols(),
    K = V_old.cols();
  W_old = clippingMX(W_old, -beta_max, beta_max);
  V_old = clippingMX(V_old, -beta_max, beta_max);
  Eigen::MatrixXd Z1(n, q + 1), ones = MatrixXd::Ones(n,1),
    oneK = MatrixXd::Ones(K,K), onepK= MatrixXd::Ones(P, K),
    W_new = W_old, V_new = V_old;
  Z1 << ones, Z;
  Eigen::MatrixXd pi_old = calcProb(Z1, V_old);
  Eigen::MatrixXd r_i= Eigen::MatrixXd::Ones(n*L, K),
    r_i_s = Eigen::MatrixXd::Ones(n*L, K);
  NumericMatrix lambda10, r_i1, r_i_s1, r_i_temp, lambda1_temp;
  NumericVector lambda2_fac = lambda2/double(L);
  Eigen::MatrixXd lambda11, lambda1X, NEMoE_step_res(P+L+q+1,K),
  NEMoE_step_res1(P+L+q+1+n*L + P, K);
  Eigen::MatrixXd pen_fac(P,K);
  lambda1X = as<Eigen::Map<Eigen::MatrixXd>>(lambda1);

  NumericVector start_X(L), end_X(L), start_W(L), end_W(L);
  for(int l=0; l<L;l++){
    if(l == 0){
      start_X[l] = 0;
      start_W[l] = 0;
      end_X[l] = seg[l];
      end_W[l] = seg[l] + 1;
    }else{
      start_X[l] = end_X[l - 1];
      start_W[l] = end_W[l - 1];
      end_X[l]  = start_X[l] + seg[l];
      end_W[l]  = start_W[l] + seg[l] + 1;
    }
  }

  for(int l = 0; l<L; l++){
    r_i.block(l * n, 0, n, K) =
      compLikeli(X.block(0, start_X[l], n, seg[l]),
                  Z, y.column(l),
                  W_old.block(start_W[l], 0, seg[l] + 1, K), V_old);

  }
  r_i = clippingMX(r_i, 1e-12, 1);
  r_i_s = r_i*oneK;
  r_i = ((r_i.array())/(r_i_s.array())).matrix();
  r_i1 = wrap(r_i);
  r_i_s1 = wrap(r_i_s);
//  pen_fac = calcPen_fac1(X, seg, Z, y, W_old, V_old, adapt);
  lambda10 = wrap(lambda1X.array() * pen_fac0.array());

  NEMoE_step_res = NEMoE_EM1(X, seg, Z, y, r_i1, lambda10,
                             lambda2_fac, alpha1, alpha2, beta_max,
                             EM_opt, gamma);
  W_new = NEMoE_step_res.block(0, 0, P+L, K);
  V_new = NEMoE_step_res.block(P+L, 0, q+1, K);

  pen_fac = calcPen_fac1(X, seg, Z, y, W_new, V_new, adapt);

  if(btr){

    for(int l = 0; l < L; l++){
      r_i_temp = wrap(r_i.block(l * n, 0, n, K));
      lambda1_temp = wrap(lambda1X.block(start_X[l],0, seg[l],K));
      W_new.block(start_W[l], 0, seg[l] + 1, K) = NEMoE_btrW(
        X.block(0, start_X[l], n, seg[l]), Z, y.column(l),
        r_i_temp, lambda1_temp,
        alpha1.row(l), W_new.block(start_W[l], 0, seg[l] + 1, K),
        W_old.block(start_W[l], 0, seg[l] + 1, K), V_new,
        pen_fac.block(start_X[l],0, seg[l],K),
        pen_fac0.block(start_X[l],0, seg[l],K), adapt);
    }

    V_new = NEMoE_btrV(X, Z, y, seg, r_i1, lambda1, alpha1,
                       lambda2, alpha2, W_new, V_new, V_old, adapt, PLL0);
  }

  if(adapt){
    pen_fac = calcPen_fac1(X, seg, Z, y, W_new, V_new, adapt);
  }
  NEMoE_step_res1 << W_new, V_new, r_i, pen_fac;

  return NEMoE_step_res1;
}

//' Cpp function: fitting NEMoE paramters.
//' (observed likelihood, penalized likelihood,
//' complete likelihood, penalized complete likelihood)
//' in all levels.
//' @param X an aggregated data matrix (n*P) of input in experts network
//' (rbind of all levels input).
//' @param Z a data matrix (n*q) of input in gating network.
//' @param y a vector of response (n).
//' @param K A number of latent class.
//' @param seg an integer vector of length in each level (L).
//' @param W_init aggregated initial parameters in experts network
//' (cbind of all levels) ((P+L)*K).
//' @param V_init initial parameters in gating network((q+1)*K).
//' @param lambda1 aggregated penalty lambda parameters in experts network(P*K).
//' @param lambda2 penalty parameters lambda in gating network(q*1).
//' @param alpha1 penalty parameters alpha in experts network(L*K).
//' @param alpha2 penalty parameter alpha in gating network(1).
//' @param beta_max A number of maximal of coefficients to avoid divergence of
//' during the fitting.
//' @param EM_opt A integer indicate methods for EM algorithm.
//'  0="EM", 1 = "CEM", 2 = "SEM", 3= "SAEM", 4="GEM".
//' @param itmax maximal numbers of iteration in fitting NEMoE.
//' @param itmin minimal numbers of iteration in fitting NEMoE.
//' @param adapt A boolean variable indicates whether to use adaptive NEMoE.
//' @param btr A boolean variable indicates whether to use backtracking during
//' fitting NEMoE.
//' @param stop_all A boolean variable indicates whether to stop by
//' (likelihood converge)&(parameters converge)
//' @param verbose A boolean variable indicates whether to show PLL
//'  during each iteration.
//' @param early_stop A boolean variable indicates whether to stop when one
//' latent class have select zero variables (to save time).
//' @return A matrix of fitting result.
//'
// [[Rcpp::export]]
Eigen::MatrixXd fitNEMoE0(Eigen::MatrixXd X, NumericVector seg,
                          Eigen::MatrixXd Z, NumericMatrix y, int K,
                          NumericMatrix lambda1, NumericVector lambda2,
                          NumericMatrix alpha1, double alpha2,
                          Eigen::MatrixXd V_init, Eigen::MatrixXd W_init,
                          double beta_max, int EM_opt, int itmax, int itmin,
                          bool adapt, bool btr, bool stop_all, bool verbose,
                          bool early_stop = false){

  int n = Z.rows(), P = X.cols(), q = Z.cols(), L = seg.size();
  Eigen::MatrixXd NEMoE_step_res, lambda1X(P,K), lambda1_eff1(P,K),
  pen_fac(P,K), W(P+L,K), V(q+1,K), r_i(n*L,K), r_i0(n*L,K), W0 = W_init,
  V0 = V_init, onepK= MatrixXd::Ones(P, K);
  NumericVector PLL;
  double PLL0, PLL1, gamma;
  bool cond_itmax, cond_itmin, cond_V, cond_W, cond_loss, flag = true, stop_cond = false;
  int it = 1, nz;
  Eigen::MatrixXd pen_fac0 = calcPen_fac1(X, seg, Z, y, W_init, V_init, adapt);
  NumericMatrix lambda1_eff;
  lambda1X = as<Eigen::Map<Eigen::MatrixXd>>(lambda1);
  lambda1_eff1 = (lambda1X.array() * pen_fac0.array()).matrix();
  lambda1_eff = wrap(lambda1_eff1);

  PLL = PLL_obs(X, Z, y, seg, W0, V0, lambda1_eff, lambda2, alpha1, alpha2, true);
  PLL0 = PLL[L];

  NumericVector start_X(L), end_X(L), start_W(L), end_W(L);
  for(int l=0; l<L;l++){
    if(l == 0){
      start_X[l] = 0;
      start_W[l] = 0;
      end_X[l] = seg[l];
      end_W[l] = seg[l] + 1;
    }else{
      start_X[l] = end_X[l - 1];
      start_W[l] = end_W[l - 1];
      end_X[l]  = start_X[l] + seg[l];
      end_W[l]  = start_W[l] + seg[l] + 1;
    }
  }

  while(flag){
    gamma = double(1)/double(it);
    try{
      NEMoE_step_res = NEMoE_step(X, seg, Z, y, lambda1, lambda2, alpha1, alpha2,
                                  V0, W0, beta_max, EM_opt, pen_fac0, adapt,
                                  gamma, btr, PLL0);
    }catch(...){
      flag = false;
    }

    W = NEMoE_step_res.block(0, 0, P+L, K);
    V = NEMoE_step_res.block(P+L, 0, q+1, K);
    r_i = NEMoE_step_res.block(P+L+q+1, 0, n*L,K);
    pen_fac = NEMoE_step_res.block(P+L+q+1+n*L, 0, P, K);
    lambda1_eff1 = (lambda1X.array() * pen_fac.array()).matrix();
    lambda1_eff = wrap(lambda1_eff1);
    PLL = PLL_obs(X, Z, y, seg, W, V, lambda1_eff, lambda2, alpha1, alpha2, true);
    PLL1 = PLL[L];
    cond_itmax = (it > itmax);
    cond_itmin = (it >= itmin);
    cond_V = ((((V - V0).array().square().sum())/(V0.array().square().sum())) < 1e-7);
    cond_W = ((((W - W0).array().square().sum())/(W0.array().square().sum())) < 1e-7);
    cond_loss = ((PLL1 - PLL0) < 1e-6);

    if(early_stop){
      for(int l = 0; l<L; l++){
        for(int k = 0; k<K; k++){
          nz = calc_nz_row(W.block(start_W[l], 0, seg[l] + 1, K).col(k));
          if(nz == 0){
            flag = false;
            W0 = W;
            V0 = V;
            PLL0 = PLL1;
            r_i0 = r_i;
            pen_fac0 = pen_fac;
            break;
          }
        }
      }

    }


    if(cond_itmin){
      if(cond_itmax){
        flag = false;
      }else{
        if(stop_all){
          stop_cond = cond_V || cond_W || cond_loss;
          if(stop_cond){
            flag = false;
          }
        }else{
          stop_cond = cond_V && cond_W && cond_loss;
          if(stop_cond){
            flag = false;
          }
        }
      }
    }
    if(flag){
      if(verbose){
        if(it == 1){
          Rcout<< "Fitting NEMoE.... \n"<<"it: "<<0<<", PLL0:"<<PLL0 <<"\n";
          Rcout<<"it: "<<it<<", PLL:"<<PLL1 <<"\n";
        }else{
          Rcout<<"it: "<<it<<", PLL:"<<PLL1 <<"\n";
        }
      }
      W0 = W;
      V0 = V;
      PLL0 = PLL1;
      r_i0 = r_i;
      pen_fac0 = pen_fac;
      it = it + 1;
    }

  }
  NEMoE_step_res << W0, V0, r_i0, pen_fac0;
  return NEMoE_step_res;
}


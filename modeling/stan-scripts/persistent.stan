/*
Title: Model 27: persistent CRP, Gaussian likelihood, one size prior
Date: 16 Jan 2023
Authors: Dan-Mircea Mirea & Yeon Soon Shin
Description: individual Stan model for the QIP (microbes) task; each individual parameter is drawn from a prior gamma distribution with parameters estimated based on a hierarchical fit of a random subset of the data (n=200);

Parameters (3):
-CRP concentration: alpha
-stay probability: eta
-size priors: sigma (variance of Gaussian generalization functions - same for both)
*/


data {
    int<lower=0> N;
    int<lower=0> NS;
    int<lower=0> K;
    int<lower=0> trial[N];
    int<lower=0> resp[N];
    int<lower=0> nk[N,K];
    int<lower=0> kmax[N];
    int<lower=0> block_[N];
    int<lower=1> sub_nos[N];
    
    real num_sq_dist[N,K];
    real len_sq_dist[N,K];
}

parameters {
    ////sample params////
	   
    // concentration
    real<lower=0> alpha_ind;

    // persistence
    real<lower=0, upper=1> eta_ind;

    // size prior - variance
    real<lower=0> sigma_ind_pre;
}

transformed parameters {
    real sigma_ind;

    // scaling the values to the 13-level space
    sigma_ind = sigma_ind_pre*13*13; 

}

model {
    //// define ////
    
    int sub = 0;
    
    vector[K] prior_k;
    vector[K] lik_k_num;
    vector[K] lik_k_len;
    vector[K] lik_k;
    vector[K] posterior_k;
    vector[K] w = rep_vector(0.0, K);;

    ////sample params individually from group-level gamma distributions////
    
    // concentration
    alpha_ind ~ gamma(0.656,0.216);

    // persistence
    eta_ind ~ beta(3.54,2.71);

    // size prior means for both dimensions
    sigma_ind_pre ~ gamma(2.76,31.5);

    ////trial-by-trial prior,lik,posterior calculation////
    
    for (i in 1:N) {

        prior_k = rep_vector(0.0, K);
	lik_k_num = rep_vector(0.0, K);
	lik_k_len = rep_vector(0.0, K);
        lik_k = rep_vector(0.0, K);
        posterior_k = rep_vector(0.0, K);
      
        // if first trial, reload subject number and assign probabilites to 1 - NO DRAWING from posterior!
        if (trial[i] == 1) {

	  sub = sub_nos[i];
	  
          prior_k[1] = 1.0;
          lik_k[1] = 1.0;
          posterior_k[1] = 1.0;
	  
        } else {
          // if not first trial, do the regular things

	  // prior for old clusters
	  prior_k = to_vector(nk[i,:]);

	  // likelihood for old clusters
	  for (k in 1:(kmax[i])){
	      lik_k_num[k] = exp(-num_sq_dist[i,k]/(2*sigma_ind));
	      lik_k_len[k] = exp(-len_sq_dist[i,k]/(2*sigma_ind));
	  }
	  lik_k = lik_k_num .* lik_k_len;

	  // prior for new cluster is equal to concentration (unnormalised)
	  prior_k[kmax[i]+1] = alpha_ind;

	  // likelihood for new cluster: the first stimulus of a cluster automatically has likelihood 1 
	  lik_k[kmax[i]+1] = 1.0;

	  // normalize prior and add persistence
    	  prior_k = (1-eta_ind)*prior_k/sum(prior_k);
    	  prior_k[resp[i-1]] += eta_ind;

    	  // normalize likelihood 
    	  lik_k = lik_k/sum(lik_k);
    
	  // compute and normalise posterior
    	  posterior_k = prior_k .* lik_k;
    	  posterior_k = posterior_k/sum(posterior_k);

    	  // statement indicating that response is drawn from categorical distribution according to posterior
    	  resp[i] ~ categorical(posterior_k);
        
       } // end for if-else statement about first trial
   
  } // end for loop through all trials
  
}


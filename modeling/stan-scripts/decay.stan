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
    real<lower=0> lambda_ind;

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
    alpha_ind ~ gamma(0.763,5.97);

    // persistence
    lambda_ind ~ gamma(1.53,5.80);

    // size prior means for both dimensions
    sigma_ind_pre ~ gamma(2.37,26.3);

    ////trial-by-trial prior,lik,posterior calculation////
    
    for (i in 1:N) {

        prior_k = rep_vector(0.0, K);
	lik_k_num = rep_vector(0.0, K);
	lik_k_len = rep_vector(0.0, K);
        lik_k = rep_vector(0.0, K);
        posterior_k = rep_vector(0.0, K);
      
        // if first trial, reload subject number and reset probabilities
        if (trial[i] == 1) {

	  sub = sub_nos[i];

	  // reset decayed evidence counter
	  w = rep_vector(0.0, K);
	  
          prior_k[1] = 1.0;
          lik_k[1] = 1.0;
          posterior_k[1] = 1.0;
	  
        } else {
          // if not first trial, do the regular things

	  // decay recency
	  w = w*exp(-lambda_ind);

	  // prior for old clusters is equal to decayed evidence counter (unnormalised)
	  prior_k = w;

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
        
       } // end for if-else statement about first trial
    
    // normalize prior and likelihood
    prior_k = prior_k/sum(prior_k);
    lik_k = lik_k/sum(lik_k);
    
    // compute and normalise posterior
    posterior_k = prior_k .* lik_k;
    posterior_k = posterior_k/sum(posterior_k);

    // statement indicating that response is drawn from categorical distribution according to posterior
    resp[i] ~ categorical(posterior_k);
    
    // update recency in evidence counter
    w[resp[i]] += + 1;

  } // end for loop through all trials
  
}

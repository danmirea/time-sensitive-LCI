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

    // sticky boost
    real<lower=0> beta_ind;

    // size prior - variance
    real<lower=0> sigma_ind_pre;

}

transformed parameters {
    real sigma_ind;

    for (sub in 1:NS){
	sigma_ind = sigma_ind_pre*13*13;
	
    }
}

model {
    //// define ////
    int sub = 0;
    vector[K] prior_k;
    vector[K] lik_k_num;
    vector[K] lik_k_len;
    vector[K] lik_k;
    vector[K] posterior_k;

    ////sample params////
    
    // concentration
    alpha_ind ~ gamma(0.830, 0.259);

    // sticky boost 
    beta_ind ~ gamma(1.63, 0.0234);

    // size prior means for both dimensions
    sigma_ind_pre ~ gamma(2.82, 31.5);

    ////trial-by-trial prior,lik,posterior calculation////
    
    for (i in 1:N) {

        prior_k = rep_vector(0.0, K);
	lik_k_num = rep_vector(0.0, K);
	lik_k_len = rep_vector(0.0, K);
        lik_k = rep_vector(0.0, K);
        posterior_k = rep_vector(0.0, K);
      
        // if first trial, reset probabilities
        if (trial[i] == 1) {
          sub = sub_nos[i];

          prior_k[1] = 1.0;
          lik_k[1] = 1.0;
          posterior_k[1] = 1.0;
	  
        } else {
          // if not first trial, do the regular things

	  // prior for old clusters
	  prior_k = to_vector(nk[i,:]);
	  prior_k[resp[i-1]] += beta_ind;

	  // likelihood for old clusters
	  for (k in 1:(kmax[i])){
	      lik_k_num[k] = exp(-num_sq_dist[i,k]/(2*sigma_ind));
	      lik_k_len[k] = exp(-len_sq_dist[i,k]/(2*sigma_ind));
	  }
	  //print(lik_k_num);
	  //print(lik_k_len);

	  lik_k = lik_k_num .* lik_k_len;

	  // prior for new cluster
	  prior_k[kmax[i]+1] = alpha_ind;

	  // likelihood for new cluster: the first stimulus has likelihood 1 
	  lik_k[kmax[i]+1] = 1.0;
        
       } // end if: first trial
    
    
    // normalize each
  
    //print("Prior is");
    //print(prior_k);
    //print("Likelihood is");
    //print(lik_k);
    
    prior_k = prior_k/sum(prior_k);
    lik_k = lik_k/sum(lik_k);

    //print("Post-norm prior is");
    //print(prior_k);
    //print("Post-norm likelihood is");
    //print(lik_k);
    
    // posterior
    posterior_k = prior_k .* lik_k;

    //print(posterior_k);

    posterior_k = posterior_k/sum(posterior_k);
   
    //print(posterior_k);
    
    resp[i] ~ categorical(posterior_k);
    

  } // end for loop: vector row for each trial
  
}


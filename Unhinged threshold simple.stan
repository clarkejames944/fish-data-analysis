data {
    int<lower=1> N_EBSm;
    int<lower=1> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    vector[Ngroups] bp;
    real<lower=0> mu_bp;
    real <lower=0, upper=100> sigma_bp;
    
    real<lower=0, upper=100> error;    ///Leave a constant variance across fish individuals for now
    
    real interceptbefore[Ngroups]; //This should be different now as I have 2 different values within each group
    real interceptafter[Ngroups];
  

    real beta[Ngroups]; //slope will be common for each group
}

transformed parameters {
    real x2[N_EBSm];
    real x3[N_EBSm];
    vector[N_EBSm] yhat;

        for (i in 1:N_EBSm) {
        if (prev[i] < bp[fishID[i]]){
            x2[i] = 1;
         } else {
            x2[i] = 0;
         }
        }

        for (i in 1:N_EBSm) {
        if (prev[i] < bp[fishID[i]]){
            x3[i] = 0;
         } else {
            x3[i] = 1;
         }
        }
    
    
        for (i in 1:N_EBSm) {
            yhat[i] = (interceptbefore[fishID[i]] * x2[i]) + (interceptafter[fishID[i]]*x3[i]) + beta[fishID[i]] * (prev[i]-bp[fishID[i]]);
      }
}

model {
    sigma_bp ~ uniform(0, 100);
    bp ~ normal(mu_bp, sigma_bp);
    
    interceptbefore ~ normal(0, 100); //intercept varies but is kept constant for the mean
    interceptafter ~ normal(0, 100);
    
    beta ~ normal(0.5, 10); ///beta is constant 

    error ~ uniform(0, 100);    

   
    oto_size ~ normal(yhat, error);
    
}


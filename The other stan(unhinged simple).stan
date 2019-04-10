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
    
    real<lower=0> error;    ///Leave a constant variance across fish individuals for now
    
    real interceptbefore[Ngroups]; //This should be different now as I have 2 different values within each group
    real interceptafter[Ngroups];


    real beta[Ngroups]; //slope will be common for each group
}

transformed parameters {
    vector[N_EBSm] yhat;

        for (i in 1:N_EBSm) {
        if (prev[i] < bp[fishID[i]]){
            yhat[i] = interceptbefore[fishID[i]] + beta[fishID[i]] * (prev[i] - bp[fishID[i]]); //got the same slopes =beta
         } else {
            yhat[i] = interceptafter[fishID[i]] + beta[fishID[i]] * (prev[i] - bp[fishID[i]]); //intercepts should be different
         }
        }
}
model {
    mu_bp ~ normal(2, 10);
    bp ~ normal(mu_bp, 100);
    
    interceptbefore ~ normal(0, 10); //intercept varies but is kept constant for the mean
    interceptafter ~ normal(0, 10);
    
    beta ~ normal(0.5, 10); ///beta is constant 

    error ~ uniform(0, 100);    
    oto_size ~ normal(yhat, error);
}

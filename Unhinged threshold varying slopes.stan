data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    vector[Ngroups] bp;
    real mu_bp;
    real <lower=0, upper=100> sigma_bp;
    
    real<lower=0, upper=100> error;    ///Leave a constant variance across fish individuals for now
    
    vector[Ngroups] interceptbefore; //This should be different now as I have 2 different values within each group
    vector[Ngroups] interceptafter;
    real<lower=0, upper=100> sigma_int_aft;
    real<lower=0, upper=100> sigma_int_bef;

    vector[Ngroups] betabefore; //slope will be different within the groups
    vector[Ngroups] betaafter;
}

transformed parameters {
vector[N_EBSm] yhat;   //The conditional mean

    for (i in 1:N_EBSm) {
        if (prev[i] < bp[fishID[i]]){
            yhat[i] = interceptbefore[fishID[i]] + betabefore[fishID[i]] * (prev[i] - bp[fishID[i]]); //got the same slopes =beta
         } else {
            yhat[i] = interceptafter[fishID[i]] + betaafter[fishID[i]] * (prev[i] - bp[fishID[i]]); //intercepts should be different (need to do this)
         }
        }
}


model {
    bp ~ normal(mu_bp, sigma_bp);
    
    interceptbefore ~ normal(0, sigma_int_bef); //intercept varies but is kept constant for the mean
    interceptafter ~ normal(0, sigma_int_aft);
    
    betabefore ~ normal(0, 1); ///beta is constant 
    betaafter ~ normal(0,1;)
    
    error ~ uniform(0, 100);    
    oto_size ~ normal(yhat, error);
}
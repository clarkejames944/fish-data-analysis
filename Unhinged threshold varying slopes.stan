data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    vector<lower=0>[Ngroups] bp;

    real<lower=0> mu_bp;
    real<lower=0> sigma_bp;
    
    real<lower=0> error;    ///Leave a constant variance across fish individuals for now
    
    real<lower=0> sigma_int_aft;
    real<lower=0> sigma_int_bef;
    vector[Ngroups] interceptbefore; //This should be different now as I have 2 different values within each group
    vector[Ngroups] interceptafter;
 

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
    mu_bp ~ normal(0, 5);
    sigma_bp ~ cauchy(0, 10);
    bp ~ normal(mu_bp, sigma_bp);
    
    sigma_int_bef ~ cauchy(0, 10);
    sigma_int_aft ~ cauchy(0, 10);
    interceptbefore ~ normal(0, sigma_int_bef); //intercept varies but is kept constant for the mean
    interceptafter ~ normal(0, sigma_int_aft);
    
    betabefore ~ normal(0, 10); ///beta is constant 
    betaafter ~ normal(0, 10);
    
    error ~ cauchy(0, 10);

    oto_size ~ normal(yhat, error);
}
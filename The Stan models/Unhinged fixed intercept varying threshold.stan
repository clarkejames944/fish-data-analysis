data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    vector<lower=0>[Ngroups] eta; 

    real<lower=0.25, upper=4> mu_eta; //allow to find the mean value for the breakpoint
    real<lower=0, upper=100> sigma_eta; //allow to find the value of sd for the breakpoint
    
    real<lower=0, upper=100> epsilon;    
    
    real alpha1; 
    real alpha2; 
    real beta;
}

transformed parameters {
    vector[N_EBSm] yhat;
    
    for (i in 1:N_EBSm) {
        if (prev[i] < eta[fishID[i]]) {
            yhat[i] = alpha1 + beta * prev[i];
        } else {
            yhat[i] = alpha2 + beta * prev[i];
        }
    }
}

model {///make sure you have a distribution for each parameter defined 
    mu_eta ~ normal(0, 10);
    sigma_eta ~ cauchy(0, 5);
    eta ~ normal(mu_eta, sigma_eta);
    
    alpha1 ~ normal(0, 10);
    alpha2 ~ normal(0, 10);
    
    beta ~ normal(0, 10);
    
    epsilon ~ cauchy(0, 10);

    oto_size ~ normal(yhat, epsilon);
}


generated quantities {
    vector[N_EBSm] sim_oto_size;

    for (i in 1:N_EBSm){
    sim_oto_size[i] = normal_rng(yhat[i], epsilon);
    }
}
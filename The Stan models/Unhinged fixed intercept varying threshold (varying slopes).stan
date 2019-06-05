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
    real<lower=0> sigma_eta; //allow to find the value of sd for the breakpoint
    
    real<lower=0, upper=10> epsilon;    
    
    real alpha1; 
    real alpha2; 
    real beta1;
    real beta2;
}

transformed parameters {
    vector[N_EBSm] tau;
    vector[N_EBSm] yhat;

    for (i in 1:N_EBSm) {
        tau[i] = 1/(1 + exp(-100 * (prev[i] - eta[fishID[i]])));
    }
    
    for (i in 1:N_EBSm){
         yhat[i] = (alpha1 + beta1 * prev[i]) * (1-tau[i]) + (alpha2 + beta2 * prev[i]) * tau[i];
        }
}
    

model {///make sure you have a distribution for each parameter defined 
    mu_eta ~ normal(1, 0.025);
    sigma_eta ~ cauchy(0.5, 0.1);
    eta ~ normal(mu_eta, sigma_eta);
    
    alpha1 ~ normal(0.275, 0.1);
    alpha2 ~ normal(0.29, 0.1);
    
    beta1 ~ normal(1, 0.05);
    beta2 ~ normal(1, 0.05);

    epsilon ~ cauchy(0, 0.1);

    oto_size ~ normal(yhat, epsilon);
}


generated quantities {
    vector[N_EBSm] sim_oto_size;
    vector[N_EBSm] log_lik;


    for (i in 1:N_EBSm){
    sim_oto_size[i] = normal_rng(yhat[i], epsilon);
    }

    for (i in 1:N_EBSm){
        log_lik[i] = normal_lpdf(oto_size[i] | yhat[i], epsilon); 
    }
}
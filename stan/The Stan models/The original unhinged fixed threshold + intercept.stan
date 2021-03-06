data {
    int<lower=0> N_EBSm;
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0.25, upper=4> eta; 
    
    real<lower=0, upper=100> epsilon; 
    
    real alpha1; 
    real alpha2;
    real beta;
}

transformed parameters {
    vector[N_EBSm] tau; //indicator variable
    
    for (i in 1:N_EBSm) {
        if (prev[i] < eta) {
            tau[i] = 0;
        } else {
            tau[i] = 1;
        }
    }
}

model {///make sure you have a distribution for each parameter defined 
    vector[N_EBSm] yhat;
    eta ~ normal(0, 10);
    
    alpha1 ~ normal(0, 10);
    alpha2 ~ normal(0, 10);

    beta ~ normal(0, 10);
    
    epsilon ~ cauchy(0, 10);

    for (i in 1:N_EBSm){
         yhat[i] = (alpha1 + beta * prev[i]) * (1-tau[i]) + (alpha2 + beta*prev[i]) * tau[i];
        }

    oto_size ~ normal(yhat, epsilon);
}

generated quantities {
    vector[N_EBSm] sim_oto_size;

    for (i in 1:N_EBSm){
    sim_oto_size[i] = normal_rng((alpha1 + beta * prev[i]) * (1-tau[i]) + (alpha2 + beta*prev[i]) * tau[i], epsilon);
    }
}
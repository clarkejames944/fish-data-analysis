data {
    int<lower=0> N_EBSm;
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0> epsilon;
    
    real alpha;
    real beta;
}

transformed parameters {
    vector[N_EBSm] yhat;
    
    for (i in 1:N_EBSm){
         yhat[i] = alpha + beta * prev[i];
        }
}

model {
    alpha ~ normal(0, 10);
    beta ~ normal(0, 10);
    
    epsilon ~ cauchy(0, 10);

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
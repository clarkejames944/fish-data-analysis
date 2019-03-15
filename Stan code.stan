data {
    int<lower=0> N_EBSm;
    int<lower=0> Npreds;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID;
    vector[N_EBSm] prev;
    vector[N_EBSm] oto_size;
}

parameters {
    real<lower=1, upper=3> breakpoint;
    
    vector[Npreds] gamma;
    real<lower=0> sigmaint;
    real<lower=0> sigmaslope;
    real<lower=0> sigmaeps;

    vector[Ngroups] etaint; //intercept
    vector[Ngroups] etaslope1; //slope before breakpoint
    vector[Ngroups] etaslope2; // part of slope after breakpoint   
}

transformed parameters {
    vector[Ngroups] beta; //slope after breakpoint
    real kappa;

    vector[Ngroups] ranint;
    vector[Ngroups] ranslope1;
    vector[Ngroups] ranslope2;
    vector[N_EBSm] yhat;

    beta <- etaslope1 + etaslope2;
    
    for (i in 1:N_EBSm) {
        if (prev[i] < breakpoint) {
            kappa[i] = 0;
        } else {
            kappa[i] = 1;
        }
    }   

    ranint <- sigmaint * etaint;
    ranslope1 <- sigmaslope * etaslope1;
    ranslope2 <- sigmaslope * etaslope2;

    for (i in 1:N_EBSm)
        yhat[i] <- ranint[fishID[i]] + ranslope1[fishID[i]] * prev[i] + ranslope2[fishID[i]] * (prev[i] - breakpoint) * kappa[i];
}

model {
    breakpoint ~ normal(0, 1);
    etaint ~ normal(0, 1);
    etaslope1 ~ normal(0, 1);
    etaslope2 ~ normal(0, 1);
    sigmaeps ~ normal(0, 2);
    sigmaint ~ normal(0, 2);
    sigmaslope ~ normal(0, 2);

    oto_size ~ normal(yhat, sigmaeps);
}
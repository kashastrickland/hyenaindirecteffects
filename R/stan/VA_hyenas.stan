
data {
  int<lower=1> N;  // total number of observations
  int<lower=0> Y[N]; // response variable
  int<lower=1> K;  // number of predictors
  matrix[N, K] X;  // population-level design matrix

  int<lower=0> N_animal;             // number of individuals in pedigree
  int<lower=0> N_NoParents;
  int<lower=0> N_ParentsOffspring;
  int<lower=0> N_ParentsNoOffspring;
  vector[N_animal] MSV;               // Mendelian sampling variance
  int<lower=0, upper=N_animal> animal_id[N];     // animal id relating y to animal in pedigree (refers to row number of dam and sire)
  int<lower=0, upper=N_animal> dam[N_animal];       // dam id in order of animal id
  int<lower=0, upper=N_animal> sire[N_animal];      // sire id  
  int<lower=0, upper=N_animal> NoParents[N_NoParents];
  int<lower=0, upper=N_animal> ParentsOffspring[N_ParentsOffspring];
  int<lower=0, upper=N_animal> ParentsNoOffspring[N_ParentsNoOffspring];

  // data for group-level effects of mother
  int<lower=1> N_mother;  // number of grouping levels
  array[N] int<lower=1> mother_id;  // grouping indicator per observation

  // data for group-level effects of birth year
  int<lower=1> N_by;  // number of grouping levels
  array[N] int<lower=1> by_id;  // grouping indicator per observation

  // data for group-level effects of birth clan
  int<lower=1> N_clan;  // number of grouping levels
  array[N] int<lower=1> clan_id;  // grouping indicator per observation
}

transformed data{
  vector[N_animal] MSsd = sqrt(MSV);
}

parameters {
  vector[K] b;  // regression coefficients
  
  real<lower=0> shape;  // shape param for neg binomial

  real<lower=0> sd_mother;  // group-level standard deviations
  vector[N_mother] z_mother;  // standardized group-level effects
  
  real<lower=0> sd_by;  // group-level standard deviations
  vector[N_by] z_by;  // standardized group-level effects

  real<lower=0> sd_clan;  // group-level standard deviations
  vector[N_clan] z_clan;  // standardized group-level effects

  real<lower=0> sd_animal;
  vector[N_NoParents] z_A_NoParents;
  vector[N_ParentsOffspring] z_A_ParentsOffspring;
  vector[N_ParentsNoOffspring] z_A_ParentsNoOffspring;
}
transformed parameters {
  vector[N_mother] r_mother;  // actual group-level effects
  vector[N_clan] r_clan;  // actual group-level effects
  vector[N_by] r_by;  // actual group-level effects

  // prior contributions to the log posterior
  // real lprior = 0;
  
  // non centred parametrisation and adding in covariance structures 
  r_mother = sd_mother * z_mother;
  r_by = sd_by * z_by;
  r_clan = sd_clan * z_clan;
}
model {
  vector [N] mu;
  vector[N_animal] z_animal;
  real z_A_mu;
  vector[N_ParentsNoOffspring] z_A_ParentsNoOffspring_mean;
 
  int PO_i; // helpful for indexing

  // dummy parent has a BV of 0
  z_animal[1] = 0;

  // founders have BV drawn from a unit normal
  z_A_NoParents ~ normal(0, 1);
  z_animal[NoParents] = z_A_NoParents;

  // Those with parents and offspring are looped over, to find the BV of the parents and take the average. The loop ensures the parents are done before their offspring. They then have a sampling variance around the mean of their parents BVs, determined by how many parents are know, and how inbred they are which is captured in MSsd
  for (i in 1:N_ParentsOffspring)
  {
    PO_i = ParentsOffspring[i];
 
    z_A_mu = (z_animal[dam[PO_i]] + z_animal[sire[PO_i]])*0.5;
    z_A_ParentsOffspring[i] ~ normal( z_A_mu, MSsd[PO_i]);
    z_animal[PO_i] = z_A_ParentsOffspring[i];
  }
  
  // For all individual with parents but no offspring, they can be samples after the loop, but because they do not have offspring can be vectorised at end
 z_A_ParentsNoOffspring_mean = (z_animal[dam[ParentsNoOffspring]] + z_animal[sire[ParentsNoOffspring]])*0.5;
  z_A_ParentsNoOffspring ~ normal( z_A_ParentsNoOffspring_mean , MSsd[ParentsNoOffspring]);
  z_animal[ParentsNoOffspring] = z_A_ParentsNoOffspring ;

  mu = exp(
      X * b +  
      r_mother[mother_id] + 
      r_by[by_id] + 
      r_clan[clan_id] + 
     sd_animal * z_animal[animal_id]);
    
  Y ~ neg_binomial_2(mu, shape);

  z_mother ~ normal(0, 1);
  z_by ~ normal(0, 1);
  z_clan ~ normal(0, 1);
  
  b ~ normal(0, 10); 
  shape ~ student_t(3, 0, 10);
  sd_mother ~ student_t(3, 0, 10);
  sd_clan ~ student_t(3, 0, 10);
  sd_by ~ student_t(3, 0, 10);
  sd_animal ~ student_t(3, 0, 10);
}
generated quantities {
  real var_mother = sd_mother^2;
  real var_by = sd_by^2;
  real var_clan = sd_clan^2;
  real var_animal = sd_animal^2;
}

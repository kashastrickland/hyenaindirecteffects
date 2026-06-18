
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

  // data for indirect effects
  int<lower=1> N_indirect;  // number of grouping levels
  array[N_indirect] int<lower=1> indirect_id;  // individuals in column of indirect design matrix
  matrix[N,N_indirect] indirect_DM;  // design matrix linking focal individuals (rows) to other individuals (columns)

  // data for group-level effects of mother
  int<lower=1> N_mother;  // number of grouping levels
  array[N] int<lower=1> mother_id;  // grouping indicator per observation

  // data for group-level effects of birth year
  int<lower=1> N_by;  // number of grouping levels
  array[N] int<lower=1> by_id;  // grouping indicator per observation

  // data for group-level effects of birth clan
  int<lower=1> N_clan;  // number of grouping levels
  array[N] int<lower=1> clan_id;  // grouping indicator per observation
  
  // data for group-level effects of spatial overlap matrix
   int<lower=1> N_hro;  // number of grouping levels
   array[N] int<lower=1> hro_id;  // grouping indicator per observation
   matrix[N_hro, N_hro] hro;  // known covariance matrix
}

transformed data{
  vector[N_animal] MSsd = sqrt(MSV);
  matrix[N_animal,2] MSsd_mat = rep_matrix(MSsd,2);
  matrix[N_hro,N_hro] L_hro = cholesky_decompose(hro);
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

  real<lower=0> sd_hro;
  vector[N_hro] z_hro;
  
  cholesky_factor_corr[2] L_corr_A;  
  vector<lower=0>[2] sd_animal;  // group-level standard deviations
  matrix [N_NoParents,2] z_A_NoParents;  // scaled breeding value
  matrix [N_ParentsOffspring,2] z_A_ParentsOffspring;
  matrix [N_ParentsNoOffspring,2] z_A_ParentsNoOffspring;

  real<lower=0> sd_indirectE;  // group-level standard deviations
  vector[N_indirect] z_indirectE;  // standardized group-level effects
}
transformed parameters {
  vector[N_mother] r_mother;  // actual group-level effects
  vector[N_clan] r_clan;  // actual group-level effects
  vector[N_by] r_by;  // actual group-level effects
  vector[N_hro] r_hro;
  vector[N_indirect] r_indirectE;  // total indirect effects, effecting each observed individual
  // prior contributions to the log posterior
  // real lprior = 0;
  
  // non centred parametrisation and adding in covariance structures 
  r_mother = sd_mother * z_mother;
  r_by = sd_by * z_by;
  r_clan = sd_clan * z_clan;
  r_hro = sd_hro * (L_hro * z_hro);
  r_indirectE =  sd_indirectE * z_indirectE;

}
model {
  vector [N] mu;
  matrix [N_animal,2] z_animal; // full vector of scaled breeding values
  row_vector[2] z_A_mu; // holds the mean BV of parents within loop below
  matrix [N_ParentsNoOffspring,2] z_A_ParentsNoOffspring_mean; // holds the mean BV of individuals with parents but no offspring
 
  int PO_i; // helpful for indexing
  matrix[N_animal, 2] A;

  // dummy parent has a BV of 0
  z_animal[1,] = rep_row_vector(0,2);

  // founders have BV drawn from a unit normal
  to_vector(z_A_NoParents) ~ normal( 0 , 1 );
  z_animal[NoParents,] = z_A_NoParents;

  // Those with parents and offspring are looped over, to find the BV of the parents and take the average. The loop ensures the parents are done before their offspring. They then have a sampling variance around the mean of their parents BVs, determined by how many parents are know, and how inbred they are which is captured in MSsd
  for (i in 1:N_ParentsOffspring)
  {
    PO_i = ParentsOffspring[i];
 
    z_A_mu = (z_animal[dam[PO_i],] + z_animal[sire[PO_i],])*0.5;
    z_A_ParentsOffspring[i,] ~ normal( z_A_mu, MSsd[PO_i]);
    z_animal[PO_i,] = z_A_ParentsOffspring[i,];
  }
  
  // For all individual with parents but no offspring, they can be samples after the loop, but because they do not have offspring can be vectorised at end
  z_A_ParentsNoOffspring_mean = (z_animal[dam[ParentsNoOffspring],] + z_animal[sire[ParentsNoOffspring],])*0.5;
  to_vector(z_A_ParentsNoOffspring) ~ normal( to_vector(z_A_ParentsNoOffspring_mean) , to_vector(MSsd_mat[ParentsNoOffspring,]));
  z_animal[ParentsNoOffspring,] = z_A_ParentsNoOffspring ;
  
  A = z_animal * diag_pre_multiply(sd_animal, L_corr_A);

  mu = exp(X * b +  
      r_mother[mother_id] + 
      r_by[by_id] + 
      r_clan[clan_id] + 
      r_hro[hro_id] +
      A[animal_id,1] + 
      indirect_DM * (r_indirectE + A[indirect_id,2]));
    
  Y ~ neg_binomial_2(mu, shape);

  z_mother ~ normal(0, 1);
  z_by ~ normal(0, 1);
  z_clan ~ normal(0, 1);
  z_hro ~ normal(0,1);

  z_indirectE ~ normal(0, 1);
  
  b ~ normal(0, 10); 
  shape ~ student_t(3, 0, 10);
  sd_mother ~ student_t(3, 0, 10);
  sd_clan ~ student_t(3, 0, 10);
  sd_by ~ student_t(3, 0, 10);
  sd_hro ~ student_t(3, 0, 10);
  L_corr_A ~ lkj_corr_cholesky(4);
  sd_animal ~ student_t(3, 0, 10);
  sd_indirectE ~ student_t(3, 0, 10);
}
generated quantities {
  real var_mother = sd_mother^2;
  real var_by = sd_by^2;
  real var_clan = sd_clan^2;
  real var_indirectE = sd_indirectE^2;
  real var_hro = sd_hro^2;
  matrix [2,2] corr_A = multiply_lower_tri_self_transpose(L_corr_A);
  matrix [2,2] Sigma_A = quad_form_diag(corr_A, sd_animal);
}

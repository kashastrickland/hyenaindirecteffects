This is the data and code repository for the manuscript in prep _The social environment shapes fitness through indirect genetic effects in spotted hyenas_ by Kasha Strickland, Matthew Silk, Oliver Höner, Eve Davidian, Loeske Kruuk and Joel Pick

The extracted life history data derived from the full hyena database is found in the folder _data_ ("lifehistory_clans_anon.csv"), together with the pedigree ("hyenapedigree_anon.csv"), social network (SRI_temporal_5_anon.csv"), home range overlap ("hro_VI_5_anon.csv", social similarity matrices ("SRItmp_sim_jac_anon.R").

All data have been anonymised with random ID names for publication and so will not correspond directly with the hyena database.

In the folder _R_ is the main code for running all models. "IGE_models_SRI.R" calls "prepare_standata.R" directly to prep data for stan models, which are in the folder _stan_. 

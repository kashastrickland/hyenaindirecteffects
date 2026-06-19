This is the data and code repository for the manuscript in prep _Social genetic effects on fitness in wild spotted hyenas_ by Kasha Strickland, Matthew Silk, Oliver Höner, Eve Davidian, Loeske Kruuk and Joel Pick

The extracted life history data derived from the full hyena database is found in the folder _data_ ("lifehistory_clans.csv"), together with the pedigree ("hyenapedigree.csv"), social network (SRI_temporal_5.csv"), home range overlap ("hro_VI_5.csv", social similarity matrices ("SRItmp_sun,jac,csv").

In the folder _R_ is the main code for running all models. "IGE_models_SRR.R" calls "prepare_standata.R" directly to prep data for stan models, which are in the folder _stan_. 

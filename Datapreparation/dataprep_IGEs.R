
## ------------------------------------------------------------
## All code required to create data from main database
## Requires access to the database -- not provided here

## ------------------------------------------------------------
## Load packages and set directory
## ------------------------------------------------------------

require(hyenaR)
require(dplyr)
require(tidyr)
library(ggplot2)
library(lubridate)
library(data.table)
library(spatsoc)
library(asnipe)
library(igraph)
library(sf)
library(adehabitatHR)


load_package_database.full(db.path="~/Documents/hyenaR_files/database_2025_09_03.sqlite")

setwd("/Users/kstrickl/Documents/GitHub/hyenaindirecteffects")

source("./Datapreparation/wrappers.R")
source("./Datapreparation/SRI_asym.R")
source("./Datapreparation/weightedJac.R")

## ------------------------------------------------------------
## extract all raw sightings from field observations & remove NAs
## ------------------------------------------------------------

raw_obvs<-extract_database_table("sightings")
raw_obvs<-raw_obvs[!is.na(raw_obvs$date_time)&
          !is.na(raw_obvs$latitude)&
          !is.na(raw_obvs$longitude)&
          !is.na(raw_obvs$ID),]

str(raw_obvs)

raw_obvsSpat<-raw_obvs[,c("date_time","latitude","longitude","ID")]

## clean up sightings that fall outside of crater
poly_crater <- st_read("~/Library/Mobile Documents/com~apple~CloudDocs/Hyenas/group_geneticdifferences_ms/Crater shapefiles 1/Crater boundary.shp")
poly_crater <- st_transform(poly_crater,sp::CRS("+proj=longlat +datum=WGS84"))
plot(poly_crater$geometry)

## frame to sf object, add a feature indicating point intersections
crater_sf <- raw_obvsSpat %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  mutate(retained = st_intersects(geometry, poly_crater, sparse = FALSE)[,1])

crater_sf$geometry <- st_transform(crater_sf$geometry, sp::CRS("+proj=utm +zone=36 +south datum=WGS84"))

crater_sf2<-crater_sf |>
  filter(retained==TRUE) |>
  mutate(longitude = st_coordinates(geometry)[,1],
         latitude = st_coordinates(geometry)[,2]) |>
  st_drop_geometry()

## check its worked OK
plot(crater_sf2$longitude,crater_sf2$latitude)
## looks OK

str(crater_sf2)
is.data.table(crater_sf2)
crater_sfDT<-as.data.table(crater_sf2)
str(crater_sf2)

## remove sightings of individuals with very few observations
sightsH<-as.data.frame(table(crater_sfDT$ID))
sightsH5<-as.data.table(sightsH[sightsH$Freq>4,])
sightsH5<-sightsH5[,Var1:=droplevels(Var1)]

crater_sfDT5<-subset(crater_sfDT,crater_sfDT$ID%in%sightsH5$Var1)

## ------------------------------------------------------------
## greate a gbi matrix of observations
## ------------------------------------------------------------

## group observations = group in 60 minute intervals
group_times(crater_sfDT5,
            datetime = "date_time",
            threshold = "60 minutes")

## group observations = group in 50 metres
## NB spatial points should in UTM here
group_pts(crater_sfDT5,
          threshold=50,
          id="ID",
          coords=c("longitude","latitude"),
          timegroup = "timegroup")

gbi<-get_gbi(crater_sfDT5,
             group="group",
             id="ID")

## ------------------------------------------------------------
## calculate networks
## ------------------------------------------------------------

## use the GBI to calculate a HWI social network
## NB: not used for analysis
net5<-get_network(gbi,
                 data_format = "GBI",
                 association_index = "HWI")

save(net5,file="data/net5.R")
write.csv(net5,file="data/socialnetwork_5.csv")

## use the GBI to calculate an assymetric SRI network
## used to deal with different differences in lifespans

SRI<-SRI_tmp(sightings=crater_sfDT5)
diag(SRI)<-1

save(SRI,file="data/SRI_temporal5.R")
write.csv(SRI,file="data/SRI_temporal_5.csv")

## ------------------------------------------------------------
## calculate social similarity
## ------------------------------------------------------------

##weighted jaccard dissimilarity as a social similarity matrix
net5_sim.jac <- jaccardMinem(net5)
net5_sim.jac[which(is.na(net5_sim.jac))] <- 0
diag(net5_sim.jac) <- 1

##weighted jaccard dissimilarity as a social similarity matrix
SRItmp_sim.jac <- jaccardMinem(SRI)
SRItmp_sim.jac[which(is.na(SRItmp_sim.jac))] <- 0
diag(SRItmp_sim.jac) <- 1

##save
save(net5_sim.jac,file="data/net5_sim.jac.R")
write.csv(net5_sim.jac,file="data/net5_sim.jac.csv")

save(SRItmp_sim.jac,file="data/SRItmp5_sim.jac.R")
write.csv(SRItmp_sim.jac,file="data/SRItmp5_sim.jac.csv")

## ------------------------------------------------------------
## calculate UDs
## ------------------------------------------------------------

str(crater_sfDT5)
hrxydata<-SpatialPointsDataFrame(crater_sfDT5[,4:5],crater_sfDT5[,"ID"])

proj4string(hrxydata) <- CRS("+proj=utm +zone=36 ellps=WGS84")  ## for example
plot(hrxydata)##visual check to make sure all looks OK

## greate grid for kernel estimation
## with 100 x 100 grid size and a km bound around sightings in crater (for smoothing purposes)
x <- seq(min(hrxydata@coords[,1])-1000,max(hrxydata@coords[,1])+1000,by=100) # where resolution is the pixel size you desire
y <- seq(min(hrxydata@coords[,2])-1000,max(hrxydata@coords[,2])+1000,by=100)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)
plot(xy)

## estimate UDs
uds7<-kernelUD(hrxydata[,1], h = "href", grid=xy)##using individual specific kernel smoothing parameters

##extract individual specific smoothing parameters and add a 50m buffer
h_values <- sapply(uds7, function(x) x@h)
h_values<-data.frame(t(h_values))
h_values$h<-as.numeric(unlist(h_values$h)+50)
h_values$ID<-rownames(h_values)

## now re=estimate the UDs with those h-vals
optud<-list()
for (i in 1:dim(h_values)[1]){
  pnts<-hrxydata[hrxydata$ID==h_values$ID[i],]
  hopt<-h_values$h[i]
  uds_man<-kernelUD(pnts,h=hopt,grid=xy)
  optud[[i]]<-uds_man
  cat(i)
}

uddf<-unlist(optud)
class(uddf)<-"estUDm"

## ------------------------------------------------------------
## calculate home range overlap
## ------------------------------------------------------------

##g oing to use VI because its on a more interpretable scale for modelling 
hro_VI<-kerneloverlaphr(uddf,method="VI",percent = 95)##using a volume intersection overlap index

## check the overlap matrices to assess which might be most appropriate for modelling
hro_VI2<-hro_VI
corpcor::is.positive.definite(hro_VI2)
hro_VI2[hro_VI2>1]<-1 ## standardise those that share all their home range (e.g., litter mates that die young)
diag(hro_VI2)<-1.001 ## standardise the diagonal and add a little jitter to make sure its pos def
corpcor::is.positive.definite(hro_VI2)

## save UDs and hro matrix 
save(uddf,file="data/uds.R")
write.csv(hro_VI2,file="data/hro_VI_5.csv")

## ------------------------------------------------------------
## get full population pedigree
## ------------------------------------------------------------

pedH<-create_pop_pedigree()
pedH<-data.frame(pedH)

write.csv(pedH,file="data/hyenapedigree.csv")

## ------------------------------------------------------------
## get life history data
## ------------------------------------------------------------

## code adapted from To Chan et al. in prep

## LRS; maternal ID; Sex; birth_clan; birth_year; birth rank (ordinal); clan_size; inbreeding coefficient; proportion of immigrant ancestry
## for dead individuals, born into the crater after start of study period

## Life history data for all IDs 
starting_table <- create_id_starting.table() # all IDs ever observed

starting_table |>
  mutate(mother = fetch_id_id.mother.genetic(ID),
         mother_social = fetch_id_id.mother.social(ID),
         father = fetch_id_id.father(ID),
         sex = fetch_id_sex(ID),
         birth_clan = fetch_id_clan.birth(ID),
         birth_date = fetch_id_date.birth(ID),
         birth_year = year(birth_date),
         mother_rank_ord = fetch_id_rank.native(ID = mother, at = birth_date), # nb of native in clan is bottom, 1 is top
         lifespan = fetch_id_duration.lifespan(ID,unit="year"),
         left_censored = fetch_id_is.censored.left(ID),
         right_censored = fetch_id_is.censored.right(ID),
         Fped = fetch_id_inbreeding(ID),
         total_number_cubs = fetch_id_number.offspring(ID)) -> complete_table_1

## Add clan size ----------------------------------------------------------
fetch_clan_number(clan = "A", at = "2000-01-01") ## Blank run for caching purposes
complete_table_1 |>
    mutate(clan_size = fetch_clan_number(clan = birth_clan, at = birth_date)) -> complete_table_2

# Selection criteria ------------------------------------------------------
# Born inside the crater
# Not adopted
# Paternity was checked by microsats
# Could have been followed for at least 10 years
# Is not left censored 
# Is not right censored

# Note about Crater born definition: by design in the data, left-censored
# individuals are considered to be born in their clan _only_ if they were females
# or males under 4 years at the date of first observation.

complete_table_2 |> 
  mutate(## flag individuals born inside the main clans with TRUE 
    keep_mainclan_born = birth_clan %in% find_clan_name.all(main.clans = TRUE),
    ## flag non-adopted IDs with TRUE 
    keep_non_adopted = mother == mother_social, 
    ## flag paternity validated with microsats with TRUE
    keep_validated_paternity = fetch_id_paternity.checked(ID),
    ## flag IDs which could have been followed for at least 10 years with TRUE (born on or before "2015-09-01") 
    keep_followed_atleast10yrs = birth_date <= (find_pop_date.observation.last() - years(10)), # 89.50% of pop. has full lifespan
    ## flag IDs that have been born after the first day of observation (born before "1996-04-12") with TRUE
    keep_left_uncensored = !left_censored,
    ## flag IDs that have been alive on the last day of observation with TRUE
    keep_right_uncensored = !right_censored) |> 
  rowwise() |> 
  ungroup() -> complete_table
  
  
## Common data for adults
lifehistory_clanIDs <- complete_table |>
  filter(keep_mainclan_born & 
           keep_followed_atleast10yrs & 
           keep_non_adopted & 
           keep_left_uncensored &
           keep_right_uncensored) |> 
  mutate(LRS = total_number_cubs,
         across(c(sex, birth_clan, mother, birth_year), as.factor)) |> 
  select(ID, LRS, lifespan, Fped,
         sex, mother_rank_ord, clan_size, birth_clan, mother, birth_year) |> 
  drop_na(ID, sex, mother_rank_ord, clan_size, birth_clan, mother, birth_year)

nrow(lifehistory_clanIDs) # 1458

## ------------------------------------------------------------
## calculate ancestry to immigrants
## ------------------------------------------------------------

## some formatting to replace NA's in parentage with clan of birth
pedHc<-pedH
pedHc$birth_clan<-fetch_id_clan.birth(pedHc$ID)
table(pedHc$birth_clan)
pedHc[pedHc$birth_clan%in%c("B","C","R","U"),]$birth_clan<-"X"

naSire <- which(is.na(pedHc[, 2]))
naDam  <- which(is.na(pedHc[, 3]))

pedHc$damOriginal <- pedHc$dam
pedHc$sireOriginal <- pedHc$sire

## replace NA in parentage with clan
pedHc[naSire, 2] <- pedHc[naSire, "birth_clan"]
pedHc[naDam, 3] <- pedHc[naDam, "birth_clan"]

## calculate ancestry proportions
Qc <- nadiv::ggcontrib(pedHc[, 1:3], ggroups = unique(pedHc$birth_clan))
pedHc <- cbind(pedHc, Qc)

## join to main dataset
lifehistory_clanIDs2<-left_join(lifehistory_clanIDs,pedHc[c("ID","X")])

## save dataset
write.csv(lifehistory_clanIDs2,file="data/lifehistory_clans.csv",row.names = FALSE)


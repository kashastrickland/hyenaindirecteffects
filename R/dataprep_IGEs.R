
require(hyenaR)
require(dplyr)
require(tidyr)
library(ggplot2)
require(MCMCglmm)
require(brms)##either these packages should be able to run what I am after
require(ggplot2)
library(lubridate)
library(data.table)
library(spatsoc)
library(asnipe)
library(igraph)
library(sf)
library(adehabitatHR)

#now download the data - choose where to store it
#download_package_csv(csv.output.folder = "~/Documents/hyenaR_files/")
##and build the database - name it something informative with date so to track version - and provide local storage location
#build_package_database.full(db.name="hyena270326",input.folder ="~/Documents/hyenaR_files/")

##and now load database. 
load_package_database.full(db.path="~/Documents/hyenaR_files/database_2025_09_03.sqlite")
### lets make some social networks!
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Hyenas/hyena_socselection")

source("./R/wrappers.R")
source("./R/SRI_asym.R")
source("./R/weightedJac.R")

create_id_starting.table(from = find_pop_date.observation.first(),        # set new margins based on previous insights
                         to = find_pop_date.observation.last()) %>%
  mutate(clan_birth=fetch_id_clan.birth(ID),
         sex=fetch_id_sex(ID),
         DOB=fetch_id_date.birth(ID),
         DOD=fetch_id_date.death(ID),
         clan_last=fetch_id_clan.current(ID,DOD),
         #rank_birth=fetch_id_rank.birth(ID),
         LRS=fetch_id_number.offspring(ID),
         mother=fetch_id_id.mother.genetic(ID)) -> start_IDs

#########
##if using all sightings
raw_obvs<-extract_database_table("sightings")
raw_obvs<-raw_obvs[!is.na(raw_obvs$date_time)&
          !is.na(raw_obvs$latitude)&
          !is.na(raw_obvs$longitude)&
          !is.na(raw_obvs$ID),]

str(raw_obvs)

raw_obvsSpat<-raw_obvs[,c("date_time","latitude","longitude","ID")]

##clean sightings that fall outside of crater
poly_crater <- st_read("~/Library/Mobile Documents/com~apple~CloudDocs/Hyenas/group_geneticdifferences_ms/Crater shapefiles 1/Crater boundary.shp")
poly_crater <- st_transform(poly_crater,sp::CRS("+proj=longlat +datum=WGS84"))
plot(poly_crater$geometry)

# frame to sf object, add a feature indicating point intersections
crater_sf <- raw_obvsSpat %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  mutate(retained = st_intersects(geometry, poly_crater, sparse = FALSE)[,1])

crater_sf$geometry <- st_transform(crater_sf$geometry, sp::CRS("+proj=utm +zone=36 +south datum=WGS84"))

crater_sf2<-crater_sf |>
  filter(retained==TRUE) |>
  mutate(longitude = st_coordinates(geometry)[,1],
         latitude = st_coordinates(geometry)[,2]) |>
  st_drop_geometry()

##check its worked OK
plot(crater_sf2$longitude,crater_sf2$latitude)
##looks OK

##group everything per day/time? -- what we want is groups to be defined if they are seen together at the same time
str(crater_sf2)
is.data.table(crater_sf2)
crater_sfDT<-as.data.table(crater_sf2)
str(crater_sf2)

##may possibly want to remove sightings of individuals with very few observations
sightsH<-as.data.frame(table(crater_sfDT$ID))
sightsH10<-as.data.table(sightsH[sightsH$Freq>4,])
sightsH10<-sightsH10[,Var1:=droplevels(Var1)]

crater_sfDT10<-subset(crater_sfDT,crater_sfDT$ID%in%sightsH10$Var1)
sightsH
start_IDs$lifespan<-as.numeric(start_IDs$DOD-start_IDs$DOB)
colnames(sightsH)[1]<-"ID"
s_ID<-merge(sightsH,start_IDs,by="ID")
plot(s_ID$Freq,s_ID$lifespan)

#load("HyenaNetworks.RData")
group_times(crater_sfDT10,
            datetime = "date_time",
            threshold = "60 minutes")

group_pts(crater_sfDT10,
          threshold=50,
          id="ID",
          coords=c("longitude","latitude"),
          timegroup = "timegroup")

gbi<-get_gbi(crater_sfDT10,
             group="group",
             id="ID")

net10<-get_network(gbi,
                 data_format = "GBI",
                 association_index = "HWI")

netGraph_in <- graph.adjacency(net10, mode="undirected", diag=FALSE, weighted=TRUE)

save(net10,file="data/net5.R")
write.csv(net10,file="data/socialnetwork_5.csv")
net10<-as.matrix(read.csv(file="data/socialnetwork_10.csv",row.names = 1))

## -- asymmetric networks to deal with different differences in lifespans

real_ASI<-AI_asym(sightings=crater_sfDT10)
diag(real_ASI)<-1

SRI<-SRI_tmp(sightings=crater_sfDT10)

save(real_ASI,file="data/real_ASI.R")
write.csv(real_ASI,file="data/real_ASI_5.csv")
save(SRI,file="data/SRI_temporal10.R")
write.csv(SRI,file="data/SRI_temporal_10.csv")

is.positive.definite(net10)
diag(net10)<-1.01
is.positive.definite(net10)

##weighted jaccard dissimilarity as a social similarity matrix
net5_sim.jac <- jaccardMinem(net10)
net5_sim.jac[which(is.na(net5_sim.jac))] <- 0
diag(net5_sim.jac) <- 1

##for the ASI
real_ASI
corpcor::is.positive.definite(real_ASI)
diag(real_ASI)<-1
corpcor::is.positive.definite(real_ASI)

##weighted jaccard dissimilarity as a social similarity matrix
ASI_sim.jac <- jaccardMinem(real_ASI)
ASI_sim.jac[which(is.na(ASI_sim.jac))] <- 0
diag(ASI_sim.jac) <- 1

##temporal SRI
diag(SRI)<-1

##weighted jaccard dissimilarity as a social similarity matrix
SRItmp_sim.jac <- jaccardMinem(SRI)
SRItmp_sim.jac[which(is.na(SRItmp_sim.jac))] <- 0
diag(SRItmp_sim.jac) <- 1

###compare and save
plot(as.vector(net10),as.vector(net5_sim.jac))
is.positive.definite(net5_sim.jac)

##save
save(net5_sim.jac,file="data/net5_sim.jac.R")
write.csv(net5_sim.jac,file="data/net5_sim.jac.csv")

save(ASI_sim.jac,file="data/ASI_sim.jac.R")
write.csv(ASI_sim.jac,file="data/ASI_sim.jac.csv")

save(SRItmp_sim.jac,file="data/SRItmp10_sim.jac.R")
write.csv(SRItmp_sim.jac,file="data/SRItmp10_sim.jac.csv")

# TODO
## -- fit maternal genetic effect
## -- dilution effect across lifespan and within groups:
## -- discussion, or does the weighting of a network deal with this
## -- sampling from model to check the accuracy of parameters

##spatio-temporal overlap
##make a home range overlap matrix then make everyone who is not alive at the same time a 0

str(crater_sfDT10)
hrxydata<-SpatialPointsDataFrame(crater_sfDT10[,4:5],crater_sfDT10[,"ID"])

proj4string(hrxydata) <- CRS("+proj=utm +zone=36 ellps=WGS84")  ## for example
#hrxydataUTM <- spTransform(hrxydata, CRS())
plot(hrxydata)##visual check to make sure all looks OK

##greate grid for kernel estimation
##with 50 x 50 grid size and a km bound around sightings in crater (for smoothing purposes)
x <- seq(min(hrxydata@coords[,1])-1000,max(hrxydata@coords[,1])+1000,by=100) # where resolution is the pixel size you desire
y <- seq(min(hrxydata@coords[,2])-1000,max(hrxydata@coords[,2])+1000,by=100)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)
plot(xy)

uds7<-kernelUD(hrxydata[,1], h = "href", grid=xy)##using individual specific kernel smoothing parameters
h_values <- sapply(uds7, function(x) x@h)
h_values<-data.frame(t(h_values))
h_values$h<-as.numeric(unlist(h_values$h)+50)
h_values$ID<-rownames(h_values)

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
plot(uddf$`L-366`)
#uds_buffered<-kernelUD(hrxydata[,1], h = h_values$h, grid=xy)##using individual specific kernel smoothing parameters

##going to use VI because its on a more interpretable scale for modelling 
hro_VI<-kerneloverlaphr(uddf,method="VI",percent = 95)##using a volume intersection overlap index
#hro_UDOI<-kerneloverlaphr(uds7,method="UDOI",percent = 95)##estimates a scaled probability of using same space

##check the overlap matrices to assess which might be most appropriate for modelling
hro_VI2<-hro_VI
diag(hro_VI2)<-1.001##standardise the diagonal and add a little to make sure its pos def
corpcor::is.positive.definite(hro_VI2)
hro_VI2[hro_VI2>1]<-1##standardise those that share all their home range (i.e. litter mates that die young)
diag(hro_VI2)<-1.001##standardise the diagonal and add a little to make sure its pos def
corpcor::is.positive.definite(hro_VI2)

#hro_UDOI2<-hro_UDOI
#diag(hro_UDOI2)<-1.01##standardise the diagonal and add a little to make sure its pos def
#corpcor::is.positive.definite(hro_UDOI2)
mat_check<-data.frame(hro=gdata::unmatrix(hro_VI),
                      net=gdata::unmatrix(net10))

##save hro matrix for models later
save(uddf,file="data/uds.R")
save(hro_VI2,file="data/hro_VI2.R")
write.csv(hro_VI2,file="data/hro_VI_5.csv")

#save(hro_UDOI,file="hro_UDOI.R")

#list.files()
#load("hro.R")
#load("uds7.R")

##get full population pedigree
pedH<-create_pop_pedigree()
pedH<-data.frame(pedH)

write.csv(pedH,file="data/hyenapedigree.csv")

##get life history data
##LRS; maternal ID; Sex; birth_clan; birth_year; borth rank (ordinal); clan_size; inbreeding coefficient; proportion of immigrant ancestry
##for dead individuals, born into the crater after start of study period

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

## Here we add 'keep columns' to make the filtering processes clearer.
## We will filter the data right before modeling/plotting, etc.

# 'TRUE' in these columns means: to keep!

## Keep column explanation ------------------------------------------------

# We have 'keep columns' for the following filtering: 
#  1. Born inside the crater
#  2. Not adopted
#  5. Paternity was checked by microsats
#  7. Could have been followed for at least 10 years (our choice for the manuscript, see extract_pop_data_companion_Rscripts/check_censorship.R)
# 12. Is not left censored (our choice for the manuscript)
# 13. Is not right censored

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

##some formatting to replace NA's in parentage with clan of birth
pedHc<-pedH
pedHc$birth_clan<-fetch_id_clan.birth(pedHc$ID)
table(pedHc$birth_clan)
pedHc[pedHc$birth_clan%in%c("B","C","R","U"),]$birth_clan<-"X"

naSire <- which(is.na(pedHc[, 2]))
naDam  <- which(is.na(pedHc[, 3]))

pedHc$damOriginal <- pedHc$dam
pedHc$sireOriginal <- pedHc$sire

##replace NA in parentage with clan
pedHc[naSire, 2] <- pedHc[naSire, "birth_clan"]
pedHc[naDam, 3] <- pedHc[naDam, "birth_clan"]

##
Qc <- nadiv::ggcontrib(pedHc[, 1:3], ggroups = unique(pedHc$birth_clan))
pedHc <- cbind(pedHc, Qc)

lifehistory_clanIDs2<-left_join(lifehistory_clanIDs,pedHc[c("ID","X")])

###
write.csv(lifehistory_clanIDs2,file="data/lifehistory_clans.csv",row.names = FALSE)



##make a temporal overlap matrix
str(clanN)
clanN$DOD<-fetch_id_date.death(clanN$id)
clanN$DOB<-fetch_id_date.birth(clanN$id)
clanN$DOD[is.na(clanN$DOD)]<-max(as.Date(raw_obvs$date_time))

tempOLv<-list()
for(d in 1:nrow(clanN)){
  tmpList<-list()
  for (f in 1:nrow(clanN)){
    if(!is.na((min(clanN$DOD[d],clanN$DOD[f]))-(max(clanN$DOB[d],clanN$DOB[f])))) {
      td<-as.integer((min(clanN$DOD[d],clanN$DOD[f]))-(max(clanN$DOB[d],clanN$DOB[f])))
    } else if (clanN$DOD[d]<clanN$DOB[f]){
      td=0
    } else {
      td=clanN$DOD[d]-clanN$DOB[f]
      }
  tmp<-data.frame(alt=clanN$id[f],td=td)
  tmpList[[f]]<-tmp 
   }
  tdPairs<-do.call("rbind",tmpList)
  tdPairs$ID<-clanN$id[d]
  tempOLv[[d]]<-tdPairs
  cat(d)
  }

tempOLvDF<-do.call("rbind",tempOLv)
tempOLvDF$td[tempOLvDF$td<0]<-0
tempOLmat<-matrix(tempOLvDF$td,nrow = 2555,ncol=2555)
colnames(tempOLmat)<-rownames(tempOLmat)<-unique(tempOLvDF$ID)

tempOLmats<-tempOLmat[rownames(tempOLmat)%in%rownames(gms),colnames(tempOLmat)%in%colnames(gms)]

##create a shared maternal env matrix
create_id_starting.table() %>%
  mutate(DOB=fetch_id_date.birth(ID),
       mother=fetch_id_id.mother.genetic(ID)) -> idS

idS<-subset(idS,idS$ID%in%rownames(net15s))
idS$YOB<-format(as.Date(idS$DOB, format="%Y-%m-%d"),"%Y")
pairsM<-expand.grid(idS$ID,idS$ID)
pairsM2<-merge(pairsM,idS,by.x="Var1",by.y="ID")
pairsM3<-merge(pairsM2,idS,by.x="Var2",by.y="ID")
pairsM3$ME<-pairsM3$mother.x==pairsM3$mother.y

pairsM3$ME[is.na(pairsM3$ME)]<-FALSE
pairsM3$ME[pairsM3$ME=="FALSE"]<-0
pairsM3$ME[pairsM3$ME=="TRUE"]<-1
pairsM3$ME[pairsM3$mother.x==pairsM3$Var2|
           pairsM3$mother.y==pairsM3$Var2|
           pairsM3$mother.y==pairsM3$Var1|
           pairsM3$mother.x==pairsM3$Var1]<-1

pairsM3<-pairsM3[order(pairsM3$Var1),]
pairsM3<-pairsM3[order(pairsM3$Var2),]

##now make a matrix
ME<-matrix(pairsM3$ME,1142,1142)
diag(ME)<-1
rownames(ME)<-colnames(ME)<-idS$ID

##edit hro matrix to make 0 if no temporal overlap
hrosT<-ifelse(tempOLmats==0,0,hros)

##extract times to use in downstream analyses
crater_sfDT15$int_date<-as.numeric(as.Date(crater_sfDT15$date_time))
crater_sfDT15_time<-crater_sfDT15[,c("group","int_date")]
crater_sfDT15_time<-unique(crater_sfDT15_time)

crater_sfDT15_fem$int_date<-as.numeric(as.Date(crater_sfDT15_fem$date_time))
crater_sfDT15_fem_time<-crater_sfDT15_fem[,c("group","int_date")]
crater_sfDT15_fem_time<-unique(crater_sfDT15_fem_time)

crater_sfDT15_mal$int_date<-as.numeric(as.Date(crater_sfDT15_mal$date_time))
crater_sfDT15_mal_time<-crater_sfDT15_mal[,c("group","int_date")]
crater_sfDT15_mal_time<-unique(crater_sfDT15_mal_time)

##make a dataframe of enter and exit dates for all individuals
str(clanN)
clanN$intDOB<-as.numeric(as.Date(clanN$DOB))
clanN$intDOD<-as.numeric(as.Date(clanN$DOD))
identical(colnames(gbi),clanN$id)

gbi15F<-gbi[,colnames(gbi)%in%clanN[clanN$sex=="female",]$id]
clanNfem<-subset(clanN,clanN$id%in%unique(crater_sfDT15_fem$ID))
identical(colnames(gbi15F),clanNfem$id)

gbi15M<-gbi[,colnames(gbi)%in%clanN[clanN$sex=="male",]$id]
clanNmal<-subset(clanN,clanN$id%in%unique(crater_sfDT15_mal$ID))
identical(colnames(gbi15M),clanNmal$id)

##test if non-random associations
networks.perm <- network_permutation(gbi, data_format="GBI",
                                     association_matrix=net15, 
                                     association_index = "HWI",
                                     times=crater_sfDT15_time$int_date,
                                     enter_time = clanN$intDOB,
                                     exit_time = clanN$intDOD,
                                     permutations=1000)
meansR<-vector()
sdsR<-vector()
cvsR<-vector()

for (i in 1:nrow(networks.perm)){
  tea<-networks.perm[i,,]
  tea[is.nan(tea)]<-0
  meansR[[i]]<-mean(tea)
  sdsR[[i]]<-sd(tea)
  cvsR[[i]]<-meansR[[i]]/sdsR[[i]]
}

#meansR<-apply(networks.perm,1,function(x) {mean(na.omit(x))})
#sdsR<-apply(networks.perm,1,function(x) {sd(na.omit(x))})
#cvsR<-sdsR/meansR
cvO<-sd(net15)/mean(net15)
sum(cvsR>=cvO)
quantile(cvsR,c(0.05,0.95))
##remove permuted networks to save space
rm(networks.perm)

##females
## TODO male and female network permutations not working
networks.permF <- network_permutation(gbi15F, data_format="GBI",
                                      association_matrix=netF, 
                                      association_index = "HWI",
                                      times=crater_sfDT15_fem_time$int_date,
                                      enter_time = clanNfem$intDOB,
                                      exit_time = clanNfem$intDOD,
                                      permutations=10) ##FIXME

meansRF<-vector()
sdsRF<-vector()
cvsRF<-vector()

for (i in 1:length(networks.permF)){
  tef<-networks.permF[i,,]
  tef[is.nan(tef)]<-0
  meansRF[[i]]<-mean(tef)
  sdsRF[[i]]<-sd(tef)
  cvsRF[[i]]<-meansRF[[i]]/sdsRF[[i]]
  cat(i)
}

##meansRF<-apply(networks.permF,1,function(x) {mean(na.omit(x))})
##sdsFR<-apply(networks.permF,1,function(x) {sd(na.omit(x))})
##cvsFR<-sdsFR/meansRF
cvFO<-sd(net20F)/mean(net20F)
sum(cvsRF>=cvFO)
quantile(cvsRF,c(0.05,0.95))

##males
networks.permM <- network_permutation(gbi20M, data_format="GBI",
                                      association_matrix=net20M, 
                                      association_index = "HWI",
                                      times=raw_obvsSpatDT20_mal_time$int_date,
                                      enter_time = clanNmal$intDOB,
                                      exit_time = clanNmal$intDOD,
                                      permutations=1000)##FIXME
meansRM<-vector()
sdsRM<-vector()
cvsRM<-vector()

for (i in 1:length(networks.permM)){
  tem<-networks.permM[i,,]
  tem[is.nan(tem)]<-0
  meansRM[[i]]<-mean(tem)
  sdsRM[[i]]<-sd(tem)
  cvsRM[[i]]<-meansRM[[i]]/sdsRM[[i]]
  cat(i)
}

#meansRM<-apply(networks.permM,1,function(x) {mean(na.omit(x))})
#sdsMR<-apply(networks.permM,1,function(x) {sd(na.omit(x))})
#cvsMR<-sdsMR/meansRM
cvMO<-sd(net20M)/mean(net20M)
sum(cvsRM>=cvMO)
quantile(cvsRM,c(0.05,0.95))

##stability of social associations
# calculate lagged association rate for great tits
lar<-LAR(gbi,times=crater_sfDT15_time$int_date,10,max_time = 4000)
larF<-LAR(gbi20F,times=raw_obvsSpatDT20_fem_time$int_date,10,max_time = 4000)##FIXME gbi sex specific needs fixing
larM<-LAR(gbi20M,times=raw_obvsSpatDT20_mal_time$int_date,10,max_time = 4000)##FIXME gbi sex specific needs fixing

## plot the results
plot(lar, type='l', axes=FALSE, xlab="Time (days)", ylab="LAR", ylim=c(0,1))
axis(2)
axis(1, at=lar[,1], labels=c(1:nrow(lar)))

plot(larF, type='l', axes=FALSE, xlab="Time (days)", ylab="LAR", ylim=c(0,1))
axis(2)
axis(1, at=larF[,1], labels=c(1:nrow(larF)))

plot(larM, type='l', axes=FALSE, xlab="Time (days)", ylab="LAR", ylim=c(0,1))
axis(2)
axis(1, at=larM[,1], labels=c(1:nrow(larM)))


#########MRQAP
##check order of matrices are the same
identical(rownames(net15s),rownames(ME))
identical(rownames(net15s),rownames(tempOLmats))
identical(rownames(net15s),rownames(hrosT))
identical(rownames(net15s),rownames(rmats))
identical(rownames(net15s),rownames(gms))
gms<-gms[rownames(net15s),colnames(net15s),drop=FALSE]
identical(rownames(net15s),rownames(gms))

####Uni-v
#MEMRQ<-mrqap.dsp(net15s~ME,randomisations = 10)
#timeMRQ<-mrqap.dsp(net15s~tempOLmats,randomisations = 10)
#spaceMRQ<-mrqap.dsp(net15s~hrosT,randomisations = 10)
#grmMRQ<-mrqap.dsp(net15s~gms,randomisations = 10)
#pedMRQ<-mrqap.dsp(net15s~rmats,randomisations = 10)

###multi-v
full_MRQ<-mrqap.dsp(net15s~tempOLmats+hrosT+rmats+ME,randomisations = 1000)
full_gm_MRQ<-mrqap.dsp(net15s~tempOLmats+hrosT+gms+ME,randomisations = 1000)

save.image(file="SNAresults.RData")



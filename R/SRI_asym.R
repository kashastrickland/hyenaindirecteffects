
SRI_tmp<-function(sightings=sightings){
  z<-sightings
  z$date<-as.Date(z$date_time)
  animals<-sort(unique(z[,ID]))
  n<-length(animals)
  matnames<-list(animals,animals)
  currentmat<-matrix(c(rep(NA,n^2)),nrow=n,dimnames=matnames)
  currentdflist<-split(z, z[,ID],drop=TRUE)
  for (i in 1:nrow(currentmat)) {
    ego<-row.names(currentmat)[i]
    #Get the list element for each animal
    all_ego<-get(ego, currentdflist)
    for (j in 1:ncol(currentmat)) {
      alter<-colnames(currentmat)[j]
      #Get the list element for each animal
      all_alter<-get(alter, currentdflist)
      ##get alters sightings that overlap with time period of ego's sightings
      all_alter_tmp<-all_alter[all_alter$date>=min(all_ego$date)&all_alter$date<=max(all_ego$date)]
      #Take the intersection of the to see when the animal were in the same group
      set<-(intersect(all_ego[,group], all_alter_tmp[,group]))
      sample<-subset(all_ego, all_ego[,group] %in% set)##egos sightings with alter
      #Numerator for association index
      X<-nrow(sample)
      if (X>0){  
        X<-nrow(sample)
        #Ego without alter
        Ya<-nrow(all_ego[!all_ego$group%in%set])
        #Alter without ego during ego's time period
        Yb<-nrow(all_alter_tmp[!all_alter_tmp$group%in%set])
        #association coefficient
        SRI_tmp<-(X/(X+Ya+Yb))
        currentmat[i,j]<-SRI_tmp}
      else{currentmat[i,j]<-0}
    }
  }
  diag(currentmat)<-NA
  return(currentmat)
}


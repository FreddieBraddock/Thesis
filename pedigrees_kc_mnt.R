#Load library
library(httpgd)
library(kinship2)
library(data.table)

#Super useful website for kinship2
#http://cran.nexr.com/web/packages/kinship2/vignettes/pedigree.pdf

#open  visual interface click on link
hgd()

#read in pedigree data from mnt
pedigree_data<-fread("/mnt/data/fbraddock/pedigree_data_position_id.csv")


#Add relation (for twins)
relation = na.omit(pedigree_data[, 9:12])

#Creat a new ID which involved id, wgs and avail data
id2 <- paste(pedigree_data$id,  paste(ifelse(pedigree_data$avail==1, "*", ""), ifelse(pedigree_data$wgs==1, "WGS", "")), sep="\n")

#Create pedigree using pedigree()
pedAll_kc <- pedigree(
  famid=pedigree_data$ped,
  id=pedigree_data$id, 
  dadid=pedigree_data$father, 
  momid=pedigree_data$mother,
  sex=pedigree_data$sex, 
  affected=pedigree_data$affected,
  status=pedigree_data$status,
  relation=relation,
)

#Change id to the more inforative id2
pedAll_kc$id <- id2

#Make vector for family loop
family_id <- unique(pedigree_data$ped)


#loop through and plot all 
for(fam in family_id){
plot.pedigree(pedAll_kc[fam])
title(main=fam)
}

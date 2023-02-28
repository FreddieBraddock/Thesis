library("data.table")

# define the file path where the variant files are located
file_path <- "/mnt/data/fbraddock/familial_kc/snp/"

# define the list of individuals in the family and those unaffected
fam.list <- list("KC_family_1" = c("KC1", "KC2"), "KC_family_2" = c("KC3", "KC4", "KC5", "KC6"), 
                 "KC_family_3" = c("KC7", "KC8", "KC9", "KC10", "KC12"), "KC_family_6" = c("KC15", "KC16"), 
                 "KC_family_7" = c("KC17", "KC18"), "KC_family_8" = c("KC19", "KC20", "KC21"), 
                 "KC_family_20" = c("KC22", "KC23", "KC24"),  "KC_family_22" = c("KC25", "KC26", "KC27"), 
                 "KC_family_24" = c("KC28", "KC29", "KC30"), "KC_family_50" = c("KC32", "KC33", "KC34", "KC35", "KC36"), 
                 "KC_family_57" = c("KC38", "KC39"), "KC_family_62" = c("BUA40", "KC41", "BUA42", "KC43", "KC44", "KC45", "BUA46", "KC47", "KC48"), 
                 "KC_family_67" = c("KC49", "KC50"), "KC_family_73" = c("BUA51", "BUA52", "BUA53", "BUA54"), 
                 "KC_family_88" = c("BUA56", "KC55", "KC57"), "KC_family_325" = c("KC59", "KC61"), 
                 "KC_family_395" = c("BUA62", "KC63"))

unaffected <- c("KC22", "BUA52")


fam.list_excl_unaff <- lapply(fam.list, function(x) x[!(x %in% unaffected)])

# Convert to data.table
fam.dt <- data.table(
  family = rep(names(fam.list), lengths(fam.list)),
  individual = unlist(fam.list)
)

# mark all individuals as affected
fam.dt[, status := "affected"]

# mark KC22 and KC1 as unaffected
fam.dt[individual %in% unaffected, status := "unaffected"]
fam.dt

# create a list to store the results for each family
result_list <- list()

# loop through each unique family
for (family_name in unique(fam.dt$family)) {
  
  # subset the data.table to only include the current family
  family_dt <- data.table()
  for (individual in unique(fam.dt[family == family_name & status != "unaffected", individual])) {
    individual_dt <- fread(paste0(file_path, individual, ".GATK.snp.annovar.hg38_multianno.xls"))
    individual_dt$Chromosome_coordinates <- paste(individual_dt[,CHROM],":",individual_dt[,POS],individual_dt[,REF],">",individual_dt[,ALT], sep="")
    xx <- c('Chromosome_coordinates',individual)
     if (nrow(family_dt) <1) {
        family_dt <- rbind(family_dt, individual_dt, fill=TRUE)
     }# thhis is because family_dt cannot merge with individual_dt on the first loop without any column names
    family_dt <- merge(family_dt, individual_dt[,..xx], by="Chromosome_coordinates")
  }
  
  # get the chromosomal coordinates of the unaffected individual
  unaffected_dt <- data.table()
  for (individual in unique(fam.dt[family == family_name & status == "unaffected", individual])) {
    individual_dt <- fread(paste0(file_path, individual, ".GATK.snp.annovar.hg38_multianno.xls"))
    individual_dt$Chromosome_coordinates <- paste(individual_dt[,CHROM],":",individual_dt[,POS],individual_dt[,REF],">",individual_dt[,ALT], sep="")
    xx <- c('Chromosome_coordinates',individual)
      if (nrow(unaffected_dt) <1) {
        unaffected_dt <- rbind(unaffected_dt, individual_dt, fill=TRUE)
     }
    unaffected_dt <- merge(unaffected_dt, individual_dt[,..xx], by="Chromosome_coordinates")
  }
    print("Before removal unaffected")
  print(nrow(family_dt))
  # remove variants with the same chromosomal coordinates as those in the unaffected individual
  unaffected_individuals <- unique(unaffected_dt$Chromosome_coordinates)
  if (length(unaffected_individuals) > 0) {
    family_dt <- family_dt[!Chromosome_coordinates %in% unaffected_individuals]
  }
  print("After removal unaffected")
  print(nrow(family_dt))

#Change gnomad_genome_AF to non scientific format 
    as.numeric(format(family_dt[, AF], scientific = FALSE))

 #Split into het 
    DFM2 <- family_dt[family_dt[, INFO] %like% "AC=1",]
    #Filter MAF<0.001
    DFM3 <- DFM2[DFM2[, AF]< 0.001,]
    #filter exonicDFM4 <- DFM3[DFM3[Func == "exonic"]]
    DFM4 <- DFM3[DFM3[Func == "exonic"]]
    #split and Filter non-synonymous
    DFM5 <- DFM4[DFM4[ExonicFunc != "synonymous SNV"]]
    # filter for synonymous 
    DFM6 <- DFM4[DFM4[ExonicFunc == "synonymous SNV"]]
    #Filter intronic variants
    DFM7 <- DFM3[DFM3[Func == "intronic"]]
    #Filter all other rare variants (intergenic)
    DFM8 <- DFM3[DFM3[Func != "intronic"]]
    DFM8 <- DFM8[DFM8[Func != "exonic"]]
    #Split into homo 
    DFMhomo <- family_dt[family_dt[, INFO] %like% "AC=2",]
    #Filter MAF<0.005
    DFM9 <- DFMhomo[DFMhomo[, AF]< 0.05,]
    #filter exonic
    DFM10 <- DFM9[DFM9[Func == "exonic"]]
    #split and Filter non-synonymous
    DFM11 <- DFM10[DFM10[ExonicFunc != "synonymous SNV"]]
    # filter for synonymous 
    DFM12 <- DFM10[DFM10[ExonicFunc == "synonymous SNV"]]
    #Filter intronic variants
    DFM13 <- DFM9[DFM9[Func == "intronic"]]
    #Filter all other rare variants (intergenic)
    DFM14 <- DFM9[DFM9[Func != "intronic"]]
    DFM14 <- DFM14[DFM14[Func != "exonic"]]
 export.loc <- "/mnt/data/fbraddock/familial_kc/snp/results/"
    number.list <- list(1,2,3,4,5,6,7,8,9,10,11,12)
    dfm.list<- list(DFM3, DFM4, DFM5, DFM6, DFM7, DFM8, DFM9, DFM10, DFM11, DFM12, DFM13, DFM14)
    file.list<- c("het_rare", "het_exonic", "het_non-synonymous", "het_synonymous", "het_intronic", "het_intergenic&other","homo_rare", "homo_exonic", "homo_non-synonymous", "homo_synonymous", "homo_intronic", "homo_intergenic&other")
    for(n in number.list){
        print(n)
        write.csv(dfm.list[[n]], file=paste0(export.loc, family_name, "/", family_name,"_", file.list[[n]], "_snp.csv"), row.names=FALSE)

    }
}

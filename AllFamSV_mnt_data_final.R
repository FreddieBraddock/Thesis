library("data.table")

# define the file path where the variant files are located
file_path <- "/mnt/data/fbraddock/familial_kc/sv/"

# define the path where files will be saved
export.loc <- "/mnt/data/fbraddock/familial_kc/sv/results/"

# define the list of individuals in the family and those unaffected

unaffected <- c("KC22", "BUA52")


fam.list <- list("KC_family_1" = c("KC1", "KC2"), "KC_family_2" = c("KC3", "KC4", "KC5", "KC6"), 
                 "KC_family_3" = c("KC7", "KC8", "KC9", "KC10", "KC12"), "KC_family_6" = c("KC15", "KC16"), 
                 "KC_family_7" = c("KC17", "KC18"), "KC_family_8" = c("KC19", "KC20", "KC21"), 
                 "KC_family_20" = c("KC22", "KC23", "KC24"),  "KC_family_22" = c("KC25", "KC26", "KC27"), 
                 "KC_family_24" = c("KC28", "KC29", "KC30"), "KC_family_50" = c("KC32", "KC33", "KC34", "KC35", "KC36"), 
                 "KC_family_57" = c("KC38", "KC39"), "KC_family_62" = c("BUA40", "KC41", "BUA42", "KC43", "KC44", "KC45", "BUA46", "KC47", "KC48"), 
                 "KC_family_67" = c("KC49", "KC50"), "KC_family_73" = c("BUA51", "BUA52", "BUA53", "BUA54"), 
                 "KC_family_88" = c("BUA56", "KC55", "KC57"), "KC_family_325" = c("KC59", "KC61"), 
                 "KC_family_395" = c("BUA62", "KC63"))

fam.list_excl_unaff <- lapply(fam.list, function(x) x[!(x %in% unaffected)])


#list of corneal dystrophy genes
cd_genes <- c("TGFBI", "COL17A1", "KRT3", "KRT12", "TACSTD2", "CHST6", "UBIAD1", "DCN", "PIKFYVE", "KERA", "LUM", "DCN", "EPYC", "STS", "COL8A2", "ZEB1", "COL8A2", "SLC4A11", "ZNF469", "PRDM5", "CHRDL1", "TCF4", "OVOL2", "ZEB1", "GRHL2")

#list of keratoconus associated genes (PMID: 35141241, Hao et al, 2021) +(Hardcastle lab)
kc_genes <- unique(c("COL5A1", "MIR184", "LOX", "ZNF469", "VSX1", "COL4A3", "COL4A4", "COL5A1", "miR184", "LOX", "DOCK9", "IPO5", "SLC4A11", "SPARC", "STK24", "CAST", "IL1RN", "HKDC1", "IL17B", "PROB1", "SKP1", "ZNF469", "COL4A3", "VSX1", "COL4A4", "FLG", "TGFBI", "TIMP3", "SOD1", "GALNT14", "PCSK1", "PPIP5K2", "TSC1", "TUBA3D", "ADAMTS3", "BIRC2", "CD248", "COL6A2", "COL6A5", "FZD2", "LRP6", "MYOF", "PAK6", "PPP1R12A", "PPP3CC", "PTK6", "STX2", "VANGL1", "WNT1", "WNT16", "ZNF676", "ZNF765", "FNDC3B", "FOXO1", "HGF", "MPDZ-NFIB", "RAB3GAP1", "IMMP2L", "CSNK1E", "MAML2", "PNPLA2", "SMAD3", "STON2", "WNT10A", "ADAMTS9", "APEX1", "CAT", "FAS", "FASLG", "FEN1", "GPX1", "IL1A", "IL1B", "KCND3", "KIF26B", "LIG3", "MAP3K19", "MMP9", "POLG", "RAD51", "TF", "TIMP1", "TNF-α", "XRCC1", "WNT3", "WNT5A", "ATP1B1", "MRPS14", "CD46", "LRP1B", "FNDC3B", "ITGA2", "LOX", "TFAP2B", "COL12A1", "NDUFAF6", "ACTL7B", "COL5A1", "FBXW5", "EIF3A", "PIDD", "FAM76B", "GRIN2B", "GALNT6", "FOXO1", "KLF5", "TNFAIP8L3", "RORA", "SMAD3", "KIF1C", "ALDH3A1", "RAB11FIP4", "SKAP1", "COL1A1", "STK35", "NA", "COL6A1", "AIFM3"))

#Create vector for Hardcastle GWAS associated Loci of interest
loci_position = c("1q24.2", "1q25.1", "1p22.2", "2q22.1", "3q26.31", "5q11.2", "5q23.2", "6p12.3", "6q13", "8q22.1", "9p23", "9q31.3", "9q34.3", "10q21.1", "10q26.11", "11p15.5", "11q21", "12p13.1", "12q13.13", "13q14.11", "13q22.1", "15q21.2", "15q22.2", "15q22.33", "16q24.2", "17p13.2", "17p11.2", "17q21.32", "17q21.33", "20p13", "20q13.31", "21q21.3", "21q22.3", "22q11.21")



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

#Data.frame for summary variant info
number_variants_dt <- data.frame(Filtering = c("homo", "homo", "homo", "homo", "het", "het", "het", "het"), Func = c("Non-synonymous", "Synonymous", "Intronic", "Intergenic and other"))
  
# loop through each unique family
for (family_name in unique(fam.dt$family)) {
  # subset the data.table to only include the current family
  family_dt <- data.table()
  for (individual in unique(fam.dt[family == family_name & status != "unaffected", individual])) {
    individual_dt <- fread(paste0(file_path, individual, ".delly.hg38_multianno.xls"))
    individual_dt$Chromosome_coordinates <- paste(individual_dt[,Chr],"-",individual_dt[,Start],"-", individual_dt[,End], "-", individual_dt[,SVID], sep="")
    #substring gets rid of the chr at the begining of the coordinates value
    individual_dt[, Chromosome_coordinates := substring(Chromosome_coordinates, 4)] 
    names(individual_dt)[16] <- individual
    
    xx <- c('Chromosome_coordinates',individual)

     if (nrow(family_dt) <1) {
        family_dt <- rbind(family_dt, individual_dt, fill=TRUE)
     }# thhis is because family_dt cannot merge with individual_dt on the first loop without any column names
       
        family_dt <- merge(family_dt, individual_dt[,..xx], by="Chromosome_coordinates")
  }
  # get the chromosomal coordinates of the unaffected individual and modify for consistency
  unaffected_dt <- data.table()
  for (individual in unique(fam.dt[family == family_name & status == "unaffected", individual])) {
        individual_dt <- fread(paste0(file_path, individual, ".delly.hg38_multianno.xls"))
        individual_dt$Chromosome_coordinates <- paste(individual_dt[,Chr],"-",individual_dt[,Start],"-", individual_dt[,End], "-", individual_dt[,SVID], sep="")
        individual_dt[, Chromosome_coordinates := substring(Chromosome_coordinates, 4)] 
      if (nrow(unaffected_dt) < 1) {
        unaffected_dt <- rbind(unaffected_dt, individual_dt, fill=TRUE)
      }
     names(individual_dt)[16] <- individual
    xx <- c('Chromosome_coordinates',individual)


      unaffected_dt <- merge(unaffected_dt, individual_dt[,..xx], by = "Chromosome_coordinates")
  }

  # remove variants with the same chromosomal coordinates as those in the unaffected individual
  unaffected_individuals <- unique(unaffected_dt$Chromosome_coordinates)
  if (length(unaffected_individuals) > 0) {
    family_dt <- family_dt[!Chromosome_coordinates %in% unaffected_individuals]
  }

  # change format of gnomad_genome_AF and AF columns to numeric
  family_dt$gnomad_genome_AF <- suppressWarnings({as.numeric(family_dt$gnomad_genome_AF)})
  family_dt$AF <- suppressWarnings({as.numeric(family_dt$AF)})
  
  #Convert "NA" (missing values) to "0"
  family_dt$gnomad_genome_AF[is.na(family_dt$gnomad_genome_AF)] <- 0
  family_dt$AF[is.na(family_dt$AF)] <- 0

    # Key family_dt by the column "Chromosome_coordinates"
    setkey(family_dt, Chromosome_coordinates)
  
  #Filter for het
  DFM4 <- family_dt[family_dt[Func %in% "exonic"]]
  #split and Filter intronic
  DFM7 <- family_dt[family_dt[Func %in% "intronic"]]
  #Filter all other rare variants (intergenic)
  DFM8 <- family_dt[family_dt[!Func %in% "intronic"]]
  DFM8 <- DFM8[DFM8[!Func %in% "exonic"]]

  loci_position = c("1q24.2", "1q25.1", "1p22.2", "2q22.1", "3q26.31", "5q11.2", "5q23.2", "6p12.3", "6q13", "8q22.1", "9p23", "9q31.3", "9q34.3", "10q21.1", "10q26.11", "11p15.5", "11q21", "12p13.1", "12q13.13", "13q14.11", "13q22.1", "15q21.2", "15q22.2", "15q22.33", "16q24.2", "17p13.2", "17p11.2", "17q21.32", "17q21.33", "20p13", "20q13.31", "21q21.3", "21q22.3", "22q11.21")
  #Filter for GWAS locus
  DFM_GWASloci <- family_dt[cytoBand %in% loci_position]



  number.list <- list(1,2,3,4,5)

  dfm.list<- list(family_dt, DFM4, DFM7, DFM8, DFM_GWASloci)
  file.list<- c("rare", "exonic", "intronic", "intergenic&other", "GWAS_loci")

  for(n in number.list){
    print(n)
    write.csv(dfm.list[[n]], file=paste0(export.loc, family_name, "/", family_name,"_", file.list[[n]], "_sv.csv"), row.names=FALSE)

  }

}






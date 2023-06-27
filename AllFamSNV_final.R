library("data.table")

# define the file path where the variant files are located
file_path <- "/home/fbraddock/NAS/WGS/Re-anno_final/"
#"/home/fbraddock/NAS/Corneal_and_endothlial_RPKM_TPM/SCD_case_RNA.csv"
# define the path where files will be saved
export.loc <- "/home/fbraddock/NAS/WGS/Analysis/snv/results/"

# define the list of individuals in the family and those unaffected


unaffected <- c("KC22", "BUA52")

#fam.list <- list("KC_family_62" = c("BUA40", "KC41", "BUA42", "KC43", "KC44", "KC45", "BUA46", "KC47", "KC48"))

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


#vector of corneal dystrophy genes
cd_genes <- c("TGFBI", "COL17A1", "KRT3", "KRT12", "TACSTD2", "CHST6", "UBIAD1", "DCN", "PIKFYVE", "KERA", "LUM", "DCN", "EPYC", "STS", "COL8A2", "ZEB1", "COL8A2", "SLC4A11", "ZNF469", "PRDM5", "CHRDL1", "TCF4", "OVOL2", "ZEB1", "GRHL2")

#vector of keratoconus associated genes (PMID: 35141241, Hao et al, 2021) +(Hardcastle lab)
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
    individual_dt <- fread(paste0(file_path, individual, "/", individual, ".GATK.snp.annovar.hg38_multianno.xls"))
    individual_dt$Chromosome_coordinates <- paste(individual_dt[,CHROM],"-",individual_dt[,POS],"-", individual_dt[,REF],"-",individual_dt[,ALT], sep="")
    individual_dt[, Chromosome_coordinates := substring(Chromosome_coordinates, 4)] 

    xx <- c('Chromosome_coordinates',individual)
     if (nrow(family_dt) <1) {
        family_dt <- rbind(family_dt, individual_dt, fill=TRUE)
     }# thhis is because family_dt cannot merge with individual_dt on the first loop without any column names
        family_dt <- merge(family_dt, individual_dt[,..xx], by="Chromosome_coordinates")
  }

  # get the chromosomal coordinates of the unaffected individual and modify for consistency
  unaffected_dt <- data.table()
  for (individual in unique(fam.dt[family == family_name & status == "unaffected", individual])) {
        individual_dt <- fread(paste0(file_path, individual, "/", individual, ".GATK.snp.annovar.hg38_multianno.xls"))
        individual_dt$Chromosome_coordinates <- paste(individual_dt[,CHROM],"-",individual_dt[,POS], "-", individual_dt[,REF], "-", individual_dt[,ALT], sep="")
        individual_dt[, Chromosome_coordinates := substring(Chromosome_coordinates, 4)] 
        xx <- c('Chromosome_coordinates',individual)
      if (nrow(unaffected_dt) < 1) {
        unaffected_dt <- rbind(unaffected_dt, individual_dt, fill=TRUE)
      }
      unaffected_dt <- merge(unaffected_dt, individual_dt[,..xx], by="Chromosome_coordinates")
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

  
  #Split into het 
  DFM2 <- family_dt[family_dt[, INFO] %like% "AC=1",]
  #Filter MAF<0.001
  DFM3 <- DFM2[as.numeric(DFM2[, AF])< 0.001,]
  #Filter (gnomad_genome_AF) MAF<0.001
  DFM3 <- DFM3[as.numeric(DFM3[, gnomad_genome_AF])< 0.001,]
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

  

  #Filter for variants in genomic corneal dystrophy genes
  DFM_het_cd <- DFM5[GeneName %in% cd_genes]
  #Filter for variants in KC genes
  DFM_het_kc <- DFM5[GeneName %in% kc_genes]
  #Filter for variants in genomic region flagged in KC
  DFM_het_GWASloci <- DFM5[cytoBand %in% loci_position]
  
  #Split into homo 
  DFMhomo <- family_dt[family_dt[, INFO] %like% "AC=2",]
  #Filter (AF) MAF<0.05
  DFM9 <- DFMhomo[DFMhomo[, AF]< 0.05,]
  #Filter (gnomad_genome_AF) MAF<0.05
  DFM9 <- DFM9[DFM9[, gnomad_genome_AF]< 0.05,]
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

  #Merge rare exonic het and homo
  DFM15 <- rbind(DFM4, DFM10)


  #Create temporary table to store variant number data + change column name to family name
  print(family_name)
  family_number <- data.frame(
    column_name = c(nrow(DFM11), nrow(DFM12), nrow(DFM13), 
    nrow(DFM14), nrow(DFM5), nrow(DFM6), nrow(DFM7), nrow(DFM8))
  )

  setnames(family_number, old = "column_name",  new = family_name)

  number_variants_dt <- cbind(number_variants_dt, family_number)
  print(number_variants_dt)


  number.list <- list(1,2,3,4,5,6,7,8,9,10,11,12,13)

  dfm.list<- list(DFM3, DFM4, DFM5, DFM6, DFM7, DFM8, DFM9, DFM10, DFM11, DFM12, DFM13, DFM14, DFM15 )

  file.list<- c("het_rare", "het_exonic", "het_non-synonymous", "het_synonymous", "het_intronic", "het_intergenic&other", "homo_rare", "homo_exonic", "homo_non-synonymous", "homo_synonymous", "homo_intronic", "homo_intergenic&other", "het_homo_exonic")

  for(n in number.list){
    print(n)
    write.csv(dfm.list[[n]], file=paste0(export.loc, family_name, "/", family_name,"_", file.list[[n]], "_snp.csv"), row.names=FALSE)

  }


  write.csv(number_variants_dt, file=paste0(export.loc, "summary_number_variants", ".csv"), row.names=TRUE)
}






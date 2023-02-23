library(data.table)
#Import files for Family

sample.names <- c("KC1", "KC2", "KC3", "KC4", "KC5", "KC6", "KC7", "KC8", "KC9", "KC10", "KC12", "KC15", "KC16", "KC17", "KC18", "KC19", "KC20", "KC21", "KC22", "KC23", "KC24", "KC25", "KC26", "KC27", "KC28", "KC29", "KC30", "KC32", "KC33", "KC34", "KC35", "KC36","KC38", "KC39", "BUA40", "KC41", "BUA42", "KC43", "KC44", "KC45", "BUA46", "KC47", "KC48", "KC49", "KC50", "BUA51", "BUA53", "BUA54", "BUA56", "KC55", "KC57", "KC59", "KC61", "BUA62", "KC63")

family.names <- c("KC_family_1", "KC_family_2", "KC_family_3", "KC_family_6", "KC_family_7", "KC_family_8", "KC_family_20", "KC_family_22", "KC_family_24", "KC_family_50", "KC_family_57", "KC_family_62", "KC_family_67", "KC_family_73", "KC_family_88", "KC_family_325", "KC_family_395")


fam.list <- list("KC_family_1" = c("KC1", "KC2"), "KC_family_2" = c("KC3", "KC4", "KC5", "KC6"), "KC_family_3" = c("KC7", "KC8", "KC9", "KC10", "KC12"), "KC_family_6" = c("KC15", "KC16"), "KC_family_7" = c("KC17", "KC18"), "KC_family_8" = c("KC19", "KC20", "KC21"), "KC_family_20" = c("KC22", "KC23", "KC24"),  "KC_family_22" = c("KC25", "KC26", "KC27"), "KC_family_24" = c("KC28", "KC29", "KC30"), "KC_family_50" = c("KC32", "KC33", "KC34", "KC35", "KC36"), "KC_family_57" = c("KC38", "KC39"), "KC_family_62" = c("BUA40", "KC41", "BUA42", "KC43", "KC44", "KC45", "BUA46", "KC47", "KC48"), "KC_family_67" = c("KC49", "KC50"), "KC_family_73" = c("BUA51", "BUA53", "BUA54"), "KC_family_88" = c("BUA56", "KC55", "KC57"), "KC_family_325" = c("KC59", "KC61"), "KC_family_395" = c("BUA62", "KC63"))



#Naming isnt working, it overides the old file with the new, how do a make it save the old ones with unique names? 
DF <- list()
for(s in family.names){
    for(f in fam.list[[s]][]){
        
        print(s)
        print(f)
        file.loc <- paste("/mnt/data/fbraddock/familial_kc/snp/", f, ".GATK.snp.annovar.hg38_multianno.xls", sep = "")
        print(file.loc)
        DF[[f]] <- fread(file.loc)
           #Change gnomad_genome_AF to non scientific format 
        as.numeric(format(DF[[f]][, AF], scientific = FALSE))
           #Create Identifying column for shared variant analysis
        DF[[f]]$Chromosome_coordinates <- paste(DF[[f]][,CHROM],":",DF[[f]][,POS],DF[[f]][,REF],">",DF[[f]][,ALT], sep="")
    }
     DFM <- DF[[fam.list[[s]][1]]]
       
    for(f in fam.list[[s]][2:length(fam.list[[s]][])]){
        xx <- c('Chromosome_coordinates',f)
        DFM <- merge(DFM, DF[[f]][,..xx], by="Chromosome_coordinates")
        #Split into het 
        DFM2 <- DFM[DFM[, INFO] %like% "AC=1",]
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
        DFMhomo <- DFM[DFM[, INFO] %like% "AC=2",]
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
        #Export and save files

        export.loc <- "/mnt/data/fbraddock/familial_kc/snp/results/"
        number.list <- list(1,2,3,4,5,6,7,8,9,10,11,12)
        dfm.list<- list(DFM3, DFM4, DFM5, DFM6, DFM7, DFM8, DFM9, DFM10, DFM11, DFM12, DFM13, DFM14)
        file.list<- c("het_rare", "het_exonic", "het_non-synonymous", "het_synonymous", "het_intronic", "het_intergenic&other","homo_rare", "homo_exonic", "homo_non-synonymous", "homo_synonymous", "homo_intronic", "homo_intergenic&other")
        for(n in number.list){
            print(n)
            write.csv(dfm.list[[n]], file=paste0(export.loc, s, "/", s,"_", file.list[[n]], "_snp.csv"), row.names=FALSE)

        }
    }
    DF <- list()

}


#Need to remove the list DF[S], but its doesnt lieka lack of charaacters
#adding the df[s] <- Null leads to an error

DF[s] <- NULL


        
        dfm.list<- list(DFM3, DFM4, DFM5, DFM6, DFM7, DFM8, DFM9, DFM10, DFM11, DFM12, DFM13, DFM14)
        file.list<- c("het_rare", "het_exonic", "het_non-synonymous", "het_intronic", "het_intergenic&other","homo_rare", "homo_exonic", "homo_non-synonymous", "homo_intronic", "homo_intergenic&other")
            for(n in dfm.list){
                for(m in file.list){
                    write.csv(n, file=paste0(export.loc, s, "/", s,"_", m, "_snp.csv"), row.names=FALSE)
                }
            }
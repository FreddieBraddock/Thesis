library(data.table)


#Import files for Family
fam_location <- "/mnt/data/fbraddock/familial_kc/snp/results/"

#Merge all families (recursive is to go inside of the directory you are in) (it looks for all the exonic files in the shared directory)
fams <- list.files(fam_location, full.names = TRUE, recursive = TRUE, pattern = "rare_exonic_het_homo_snp.csv")

#We need to know where the variatns come from which family (function(x) location divided up, we are extractng the family name)
fam_names <- unlist(lapply(fams, function(x) {unlist(strsplit(x, "/"))[9]}))
fam_names

#lapply is  a list apply, itll apply the fread to every item of the first function (fams) and names it fam_dt. Its a list of data tables
fam_dt <- lapply(fams, fread)
names(fam_dt) <- fam_names



colnames(fam_dt[[1]])

ano <- colnames(fam_dt[[1]])[1:97]
ano <- setdiff(ano, "KC1")
ano <- c(ano, "fam_id")

#Create a empty matrix with same collumn headers (minus te KC1 KC2 etc... collumn names)
merged_dt <- data.table(matrix(ncol = length(ano), nrow = 0))
setnames(merged_dt, colnames(merged_dt), ano)
#Create loop (:= is to assign and creat a new collumn)
for (i in 1:length(fam_names)){
    dt_temp <- fam_dt[[i]]
    dt_temp[, fam_id := fam_names[i]]
    dt_temp <- dt_temp[, ..ano]
    merged_dt <- rbindlist(list(merged_dt, dt_temp), use.names=TRUE, fill=FALSE)
}

#Split up and make gene names for overlapping genes
split_gene_name <- data.frame(GeneName = unlist(strsplit(as.character(merged_dt$GeneName), ",")))

#need to alter merged_dt so that rows where there are overlapping genes, are duplicated and split up

split_merge_dt <- do.call(rbind, Map(data.frame, GeneName=lapply(strsplit(merged_dt$GeneName, ","), cbind), fam_id=merged_dt$fam_id))
split_merge_dt<- cbind(n=1:nrow(split_merge_dt), split_merge_dt)

#set split_merge_dt to a datatable rather than just a data frame
setDT(split_merge_dt) 

#Alter dt to correct format
all_genes <- unique(split_gene_name[,1])

gene_dt <- data.table(genes = all_genes)


#loop that runs

for (i in 1:length(fam_names)) {
    temp_fam_dt <- as.data.table(split_merge_dt[fam_id == fam_names[i], .N, by = GeneName])
    setnames(temp_fam_dt, "N", fam_names[i])
    gene_dt <- merge(gene_dt, temp_fam_dt, by.x = "genes", by.y = "GeneName", all.x = TRUE)

}

gene_dt[, number_of_fams := Reduce(`+`, lapply(.SD, function(x) {
    !is.na(x)})), .SDcols = fam_names ]

    
#How to sort data table
gene_dt <- gene_dt[order(number_of_fams, decreasing = TRUE)]

gene_dt <- gene_dt[number_of_fams != '1']

gene_dt

write.csv(gene_dt, file="/mnt/data/fbraddock/familial_kc/snp/results/KC_Family_all/homo_het_shared_snv_genes.csv", row.names=FALSE)


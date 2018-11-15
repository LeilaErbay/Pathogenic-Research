#AUTHOR : Leila Erbay
#PURPOSE: Create an R script for determining if a gene can be categorized as pathogenic 
#         with the use of a predictive model that has the best prediction accuracy based on the
#         selected variables used to classify pathogenicity of genes
#NOTES:   Analysis and Prediction is on GENE LEVEL

######################## Downloading/Loading Packages ###########################################
#install.packages("RMySQL")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("easycsv")
#install.packages("tidyverse")
install.packages('stringr')

# install.packages("XML")
# install.packages('RCurl')
# install.packages("rlist")
# install.packages("clustMixType")
# install.packages("PCAmixdata")
# install.packages("eclust")
# install.packages("StatMatch")




library(easycsv)  
library(RMySQL)
library(tidyverse)
library(stringr)

##############################################################################
########################## Importing Data ###################################
##############################################################################

## Choose the directory the contains all the necessary files
setwd(choose_dir())

### USE CONDITIONAL TO CHECK IF WD IS VALID

#Load File 1: gene_lengths.csv
gene_lengths <- read.csv(file.choose(), header = T, stringsAsFactors = FALSE)
names(gene_lengths) <- c("gene_symbol", "gene_length")


#Load File 2: gnomad_counts_rare_healthy_variant_counts
##simple base pair vs complex : missense vs nonsense 
# nonsense: frameshift, stop_gained, splice_acceptor, splice_donor)
numRareVariants <- read.csv(file.choose(), header = T, stringsAsFactors = FALSE)




# #Load File 3: Gene2PubMed data
# gene2pubmed <- read.delim(file.choose(), header = T, stringsAsFactors = F, quote = "", sep = "\t")
# gene2pubmed <- gene2pubmed[,2:3]


# Load File 3: Gene2PubMed data ---> Gene_Contr_Tbl
############################### MySQL QUERIES #############################

#USER MUST CREATE A DATABASE in MYSQL.
#Database name must be consistent with variables below.
#MySQL code to create DB (must be created in MYSQL):
# CREATE Database databaseName

#VARIABLES THAT USER NEEDS TO ALTER FOR RESPECTIVE INPUTS:
database = 'Capstone_db'
user = 'root'
pwd = 'Rajonrondo'
driver = MySQL()
host = 'localhost'
###DATABASE NEEDS TO EXIST ALREADY

### SETTING THE CONNECTION TO DATABASE
connection <- dbConnect(drv=driver, user = user, password = pwd, host = host, dbname = database)


###############################################################################################
############################# TASK 1: MYSQL INTERACTION #######################################
###############################################################################################

# 1) Count the number of observations of unique pubmed IDs that are
#     linked to each gene, so the output would be:
#     gene, number_of_unique_pubmed_ids_linked_to_the_gene
# 2) Count the relative influence of each gene in a publication -- e.g.
#     if there are 100 genes linked to a specific pubmed ID, you should end
#     up counting 1/100 for that gene, and the output would be:
#     gene, sum_of_relative_pubmed_contributions_linked_to_the_gene


########################### GENE2PUBMED ################################################
### CREATE TABLE FOR gene2pubmed 
dbSendQuery(connection, 'DROP TABLE IF EXISTS gene2pubmed;')

dbSendQuery(connection, 
            'CREATE TABLE  gene2pubmed (
            tax_num_id INT,
            gene_id INT,
            pubmed_id INT);')

### load in gene file into gene2pubmed table 
#NOTE: MAY NEED TO CHANGE FILE PATH WITHIN QUERY2
dbSendQuery(connection, 
            'LOAD DATA LOCAL INFILE 
            \'/Users/LeilaErbay/Desktop/PathogenicityResearch/gene2pubmed.txt\' 
            INTO TABLE gene2pubmed;')

### Deleting Unnecessary Row
dbSendQuery(connection, 
            'DELETE FROM gene2pubmed
            WHERE gene_id = 0;')

### Clean erroneous data from gene2pubmed
dbSendQuery(connection, 
            'DELETE FROM gene2pubmed
            WHERE pubmed_id IS NULL OR
            gene_id is NULL;')

######################## CREATE AGGREGATE DATA OF GENE2PUBMED : GENE_CONTR_TBL #################
### CREATE GENE_CONTR_TBL : Contribution for each gene to the respective publication id
dbSendQuery(connection, 'DROP TABLE IF EXISTS gene_contr_tbl;')


dbSendQuery(connection, 
            'CREATE TABLE gene_contr_tbl as
            SELECT gene_id, pubmed_id, ROUND(1 /num_occ, 4) as gene_contr
            FROM gene2pubmed
            INNER JOIN(																								
            SELECT pubmed_id as unq_id, COUNT(pubmed_id) as num_occ					
            FROM gene2pubmed 
            GROUP BY unq_id
            ORDER BY unq_id ) as g2
            ON g2.unq_id = pubmed_id
            ORDER BY pubmed_id ASC;')


gene_contr_tbl <- dbReadTable(connection, 'gene_contr_tbl', row.names=F)
gene_contr_tbl <- gene_contr_tbl[order(gene_contr_tbl$gene_id),]
rownames(gene_contr_tbl) <- NULL


#Load File 4:Variant_summary
variant_summary <- read.delim(file.choose(), header = T, stringsAsFactors = F, quote = "", sep = "\t")

###NOTE: variant_summary will only consist of grch37 assembly genes
variant_summary <- variant_summary[variant_summary$Assembly == "GRCh37",] 
#Retrive DF 5: total_gene_contr

########################### TOTAL_CONTR_TBL PER GENE #############################################
dbSendQuery(connection, 'DROP VIEW IF EXISTS total_contr;')
dbSendQuery(connection,
            'CREATE VIEW total_contr as
            SELECT
            gene_id, SUM(gene_contr) as total_gene_contr
            FROM
            gene_contr_tbl
            GROUP BY
            gene_id
            ORDER BY
            gene_id ASC;')

total_gene_contr <- dbReadTable(connection, 'total_contr', row.names=F)


########################### NUM_PUBS_PER_GENE ################################################
dbSendQuery(connection, 'DROP VIEW IF EXISTS num_pubs_per_gene;')
dbSendQuery(connection,
            'CREATE VIEW num_pubs_per_gene as 
            SELECT
            gene_id, Count(pubmed_id) as num_pubs
            FROM
            gene_contr_tbl
            GROUP BY
            gene_id
            Order By
            gene_id ASC;')
num_pubs_per_gene <- dbReadTable(connection, 'num_pubs_per_gene', row.names = F)
 

########################### AVG_GENE_CONTR : Sum_contr / num_pubs ################################################
dbSendQuery(connection,'DROP VIEW IF EXISTS avg_gene_contr;')
dbSendQuery(connection, 
            'CREATE VIEW avg_gene_contr as
            SELECT 
            gene_id, total_gene_contr / num_pubs as avg_contr
            FROM (
            SELECT
            total_contr.gene_id, num_pubs, total_gene_contr
            FROM
            num_pubs_per_gene
            JOIN
            total_contr
            ON
            total_contr.gene_id = num_pubs_per_gene.gene_id) as allInfo
            GROUP BY
            gene_id
            ORDER BY
            gene_id ASC;')
avg_gene_contr <- dbReadTable(connection, 'avg_gene_contr', row.names = F)


rm(database, host, pwd, user, driver, connection)
############################## FIND INTERSECTION OF GENES ###################################
# IMPORTANT DF:
#   gene_lengths
#   numRareVariants
#   variant_summary
#   gene_contr_tbl
#   total_gene_contr
#   num_pubs_per_gene
#   avg_gene_contr


### Select Unique Genes from Each file
#Character Symbols of genes
unique_genesFile1 <- unique(gene_lengths$gene_symbol)
unique_genesFile2 <- unique(numRareVariants$canonical_symbol)

#Number Values of Genes
unique_genesFile3 <- unique(gene_contr_tbl$gene_id) ## Same genes for avg_gene_contr, total_gene_contr, num_pubs_per_gene
unique_genesFile4 <- unique(variant_summary$GeneID)

unique_geneSymbols <- unique(variant_summary$GeneSymbol)

# Note:
# Unique genes in gene_lengths and numRareVariants = Gene Symbol in variant_summary
# Process:
# 1. Find intersection between gene_lengths and numRareVariants = U1
# 2. Find intersection between U1 and variant_summary$GeneSymbol= U2
# 3. Select out records from variant_summary that contains selected gene symbols
# 4. Find unique GeneIDs in selected variant_summary records = unique_geneID_postU2
# 5. Find intersection between unique_geneID_postU2 and gene_contr_tbl  = U3
# 6. apply list of genes from U3 to avg_gene_contr, total_gene_contr, and num_pubs_per_gene
# 7. finalize selected genes from all tables (unique geneIDs and geneSymbols)

u1 <- intersect(unique_genesFile1, unique_genesFile2)
u2 <- intersect(u1, unique_geneSymbols)

#Select gene_lengths and numRareVariants that intersect u2
selected_geneLengths <- gene_lengths[gene_lengths$gene_symbol %in% u2, ]
selected_numRareVariants <- numRareVariants[numRareVariants$canonical_symbol %in% u2,]
selected_variant_summary <- variant_summary[variant_summary$GeneSymbol %in% u2,]

unique_geneID_postU2 <- unique(selected_variant_summary$GeneID)
unique_geneSymbol_postU2 <- unique(selected_variant_summary$GeneSymbol)

u3 <- intersect(unique_geneID_postU2, unique_genesFile3)

uFin_genes <-  unique(selected_variant_summary$GeneID)
uFin_symbols <- unique(selected_variant_summary$GeneSymbol)
# finalize selection of geneIDs and geneSymbols -- 6140 geneIDs and related geneSymbols
# u3 contains finalized geneIds
# unique_geneSymbol_postU2  contains finalized geneSymbols

selected_geneLengths <- gene_lengths[gene_lengths$gene_symbol %in% uFin_symbols, ]
selected_numRareVariants <- numRareVariants[numRareVariants$canonical_symbol %in% uFin_symbols,]
selected_variant_summary <- variant_summary[variant_summary$GeneSymbol %in% uFin_symbols,]


selected_gene_contr <- gene_contr_tbl[gene_contr_tbl$gene_id %in% uFin_genes,]
selected_tot_gene_contr <- total_gene_contr[total_gene_contr$gene_id %in% uFin_genes, ]
selected_num_pubs <- num_pubs_per_gene[num_pubs_per_gene$gene_id %in% uFin_genes,]
selected_avg_gene_contr <- avg_gene_contr[avg_gene_contr$gene_id %in% uFin_genes,]

################## CLEAN ENVIRONMENT OF UNNECESSARY DATA FRAMES #############################
rm(avg_gene_contr, gene_contr_tbl, gene_lengths, num_pubs_per_gene, numRareVariants,
   total_gene_contr, variant_summary, u1, u2, u3, unique_geneID_postU2, unique_genesFile1,
   unique_genesFile2, unique_genesFile3, unique_genesFile4, uFin_genes, uFin_symbols,
   unique_geneSymbols, selected_geneIDs)

########################################################################################################
# Important DF that have selected data:
#   selected_geneLengths
#   selected_numRareVariants
#   selected_variant_summary
#   selected_gene_contr
#   selected_tot_gene_contr
#   selected_num_pubs
#   selected_avg_gene_contr

#### Checking intersections
intersectList <- list()
intersectList[[1]] <- setdiff(selected_geneLengths$geneID, selected_numRareVariants$canonical_symbol)
intersectList[[2]] <- setdiff(selected_geneLengths$geneID, selected_variant_summary$GeneSymbol)
intersectList[[3]] <- setdiff(selected_numRareVariants$canonical_symbol, selected_variant_summary$GeneSymbol)
intersectList[[4]] <- setdiff(selected_gene_contr$gene_id, selected_variant_summary$GeneID)
intersectList[[5]] <- setdiff(selected_avg_gene_contr$gene_id, selected_variant_summary$GeneID)
intersectList[[6]] <- setdiff(selected_tot_gene_contr$gene_id, selected_variant_summary$GeneID)
intersectList[[7]] <- setdiff(selected_num_pubs$gene_id, selected_variant_summary$GeneID)
intersectList[[8]] <- setdiff(selected_avg_gene_contr$gene_id, selected_variant_summary$GeneID)

for(i in 1:length(intersectList)){
  if(length(intersectList[[i]]) != 0) cat("Error in",i, "\n")
}


######## IF NO ERROR, CONTINUE; IF ERROR, FIX INTERSECTIONS #############################
readkey <- function()
{
  cat ("Press [enter] to continue if no error exists in the intersections")
  line <- readline()
}
readkey()

rm(i, intersectList)

#########################################################################################
########################## Clean and Set DFs ###########################################
#########################################################################################
rownames(selected_variant_summary) <- NULL
names(selected_variant_summary)[1] <- "AlleleID"
names(selected_variant_summary)[2] <- "variant_type"
names(selected_variant_summary)[3] <- "variant_name"
names(selected_variant_summary)[10] <- "RS_dbSNP"
names(selected_variant_summary)[11] <- "NSV_dbVar"

selected_variant_summary <- rowid_to_column(selected_variant_summary, "recordID")

###### NOTE: CLEAN OUT ID COLUMNS FROM selected_variant_summary

selected_variant_summary$LastEvaluated <- as.Date(selected_variant_summary$LastEvaluated, "%b %d, %Y")


# ###NOTE:::: NEED A BETTER WAY TO CONVERT ALL OTHER TO NONSENSE !!!
# for (i in 1:nrow(numRareVariants)){
#   numRareVariants[i,'canonical_csq'] <- ifelse(numRareVariants[i,'canonical_csq'] == "missense_variant", "M", "N")
# }


#Selecting the more meaningful variables from variant_summary
selected_var_summ_upd <- selected_variant_summary[,-c(4,7,11,12,13,14,17,18,19, 23,24,25,28,30,32)]

selected_var_summ_upd$NumPheno <- NA


selected_var_summ_upd$NumPheno <- lengths(regmatches(selected_var_summ_upd$PhenotypeList, gregexpr(",", selected_var_summ_upd$PhenotypeList)))+1


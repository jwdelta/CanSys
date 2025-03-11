# Load Required Library
library(optparse)
library(fgsea)

# Define the command line options
option_list = list(
    make_option(c("--gene_level_impact_scores"), type = "character", default = NULL, 
                help = "Path to the gene level impact scores file", metavar = "file"),
    make_option(c("--database_name"), type = "character", default = NULL, 
                help = "Name of the database to use", metavar = "name"),
    make_option(c("--database_resource"), type = "character", default = NULL, 
                help = "Path to the database resource file", metavar = "file"),
    make_option(c("--universe_genes"), type = "character", default = NULL, 
                help = "Path to the file containing all universe genes of humans", metavar = "file"),
    make_option(c("--nPermSimple"), type = "integer", default = NULL,
		help = "Number of permutations used in the permutation test to estimate P-values", metavar = "number"),
    make_option(c("--full_output_file_path"), type = "character", default = NULL, 
                help = "Path for the full output file", metavar = "path"),
    make_option(c("--signif_output_file_path"), type = "character", default = NULL, 
                help = "Path for the significant output file", metavar = "path")
)

# Create an OptionParser instance
opt_parser = OptionParser(option_list = option_list)
# Parse the command line options
opt = parse_args(opt_parser)

# Access the options using more descriptive variable names
gene_level_impact_scores_path <- opt$gene_level_impact_scores
database_name <- opt$database_name
database_resource_path <- opt$database_resource
universe_genes_path <- opt$universe_genes
nPermSimple <- opt$nPermSimple
full_output_file_path <- opt$full_output_file_path
signif_output_file_path <- opt$signif_output_file_path

# Function: Gene Set Enrichment Analysis
run_fgsea <- function (input_scores_path, database_name, database_resource_path, universe_genes_path) {
    # Read data
    input_scores <- read.table(input_scores_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    universe_genes <- readRDS(universe_genes_path)
    database_resource <- readRDS(database_resource_path)

    # Prepare input data
    valid_input_scores <- input_scores[input_scores[,2] != 0, , drop=FALSE]
    valid_input_matrix <- as.matrix(valid_input_scores)

    idTable <- database_resource$idTable
    collection <- database_resource$idList
    collection <- lapply(collection, as.character)
    collection <- lapply(collection, function(x) x[!is.na(x)])
    collection <- collection[sapply(collection, function(x) length(x) > 0)]

    result <- matrix(NA, ncol = 9, nrow = length(collection))
    colnames(result) <- c("N", "DE", "ES", 'NES', "GeneScoreMean", "PDS", "pval", "padj", "Affected genes")
    # N: number of genes in each GO/KEGG
    # DE: number of genes related to input gene list in each GO term
    rownames(result) <- names(collection)
    result[, "N"] <- unlist(lapply(collection, length))
    result[, "DE"] <- unlist(lapply(collection, function(x) sum((valid_input_matrix[,1] %in% x))))
    result[, "GeneScoreMean"] <- unlist(lapply(collection, function(x) mean(as.numeric(valid_input_matrix[valid_input_matrix[,1] %in% x, 2]), na.rm=TRUE)))
    result[, "Affected genes"] <- unlist(lapply(collection, function(x) {
        affected_genes <- valid_input_matrix[valid_input_matrix[,1] %in% x, 1]
        paste(affected_genes, collapse = ", ")}))
    result <- as.data.frame(result)
    result <- result[which(result[,"DE"] != 0),]
    collection <- collection[rownames(result)]
 
    valid_input_matrix <- rbind(valid_input_matrix, cbind(universe_genes[which(!(universe_genes %in% valid_input_matrix[,1]))],
	  		      rep(0, length(universe_genes[which(!(universe_genes %in% valid_input_matrix[,1]))]))))
    valid_input_matrix <- valid_input_matrix[order(as.numeric(valid_input_matrix[,2]), decreasing = FALSE),]
    ranks <- as.numeric(valid_input_matrix[,2])
    names(ranks) <- valid_input_matrix[,1]
    output <- fgseaMultilevel(pathways = collection, stats = ranks, scoreType = "pos", nPermSimple = nPermSimple)
    output <- as.data.frame(output)
    rownames(output) <- output$pathway
    output <- output[rownames(result),]
    result[,"ES"] <- output$ES
    result[,"NES"] <- output$NES
    result[,"PDS"] <- result[,"NES"]*as.numeric(result[,"GeneScoreMean"])
    result[,"pval"] <- output$pval
    result[,"padj"] <- p.adjust(result$pval, method = "BH")
    #result <- result[result$padj<0.25, , drop=FALSE]
    result <- merge(idTable, result, by.x = colnames(idTable)[1], by.y = "row.names")
    result <- result[order(result$padj),]
    result
}


# Main execution
full_gsea_result <- run_fgsea(gene_level_impact_scores_path, database_name, database_resource_path, universe_genes_path)
signif_gsea_result <- full_gsea_result[full_gsea_result$padj<0.25, , drop=FALSE]
signif_gsea_result <- signif_gsea_result[order(signif_gsea_result$padj),]

write.table(full_gsea_result, full_output_file_path, col.names = T, row.names = F, sep = "\t", quote = F)
write.table(signif_gsea_result, signif_output_file_path, col.names = T, row.names = F, sep = "\t", quote = F)



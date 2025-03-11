# Load Required Library
library(rGMMtest)

# Parse Command Line Arguments
argv <- commandArgs(TRUE)
expressionFilePath <- argv[1]
outputFilePath <- argv[2]

# Configuration Parameters
KS <- 5
performLogTransformation <- TRUE

# Read Expression Data
expressionData <- read.table(expressionFilePath, header = F, row.names = 1, check.names = F)
genes <- rownames(expressionData)

# Apply Log Transformation if Specified
expressionVector <- if (performLogTransformation) {
    log2(expressionData[,1] + 1)
} else {
    expressionData[,1]
}

# Configure and Run Gaussian Mixture Model (GMM)
custom.settings <- GMM_1D_opts
custom.settings$KS <- KS
mix_test <- runGMM(expressionVector, opts = custom.settings)

# Calculate Expression Cutoff
alpha <- mix_test$model$alpha
mu <- mix_test$model$mu
sigma <- mix_test$model$sigma
dist.plot <- generate_dist(expressionVector, alpha, mu, sigma, 1e4)
thresholds <- find_thr_by_params(alpha, mu, sigma, dist.plot)
cutoff <- min(na.omit(thresholds))

# Prepare Output Data Frame
GeneExprBinaryRep <- data.frame(Gene = genes, Expressed = ifelse(expressionVector < cutoff, 0, 1))

# Write Output to File
write.table(GeneExprBinaryRep, outputFilePath, sep = "\t", quote = F, col.names = T, row.names = F)

library(WGCNA)
options(stringsAsFactors = FALSE);

###---------------------- 
# Read data
print("Read Data")
expressionList = read.csv('expressionList_5xFAD', header = TRUE);

## Prepare and clean data
#Remove rows with less than 1 TPM
expressionList = expressionList[expressionList[,ncol(expressionList)]>1,]

datExpr0 = as.data.frame(t(expressionList[,-c(1)]));
names(datExpr0) = expressionList$gene_id;
rownames(datExpr0) = names(expressionList)[-c(1)];

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExpr0, verbose = 3);
#if not okay 
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## Clustering
sampleTree = hclust(dist(datExpr0), method = "average");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)

datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

collectGarbage();

save(datExpr, file = "data_input.RData")

###---------------------- 
## Modules construction
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

## Set Power
softPower = 15;
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
save(TOM, file = "TOM.RData")

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, dynamicColors, file = "Data-networkConstruction.RData")

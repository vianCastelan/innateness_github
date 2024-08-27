load("data/pca.RData")
pdf("../revision/principal_component_weights.pdf")
barplot(-pca$rotation[c(order(pca$rotation[,1])[1:10],order(-pca$rotation[,1])[1:10]),1], las=2, main = "top 10 PC1 weights")
barplot(-pca$rotation[c(order(pca$rotation[,2])[1:10],order(-pca$rotation[,2])[1:10]),2], las=2, main = "top 10 PC2 weights")
barplot(-pca$rotation[order(pca$rotation[,1]),1], main = "PC1 weights")
barplot(-pca$rotation[order(pca$rotation[,2]),2], main = "PC2 weights")
dev.off()



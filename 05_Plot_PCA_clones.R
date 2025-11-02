
library("SNPRelate")
vcf.fn <- "Daphnia_Helene_DmagnaLRV01.mQ100.SNP.NoMissing.vcf.tozip.gz"
snpgdsVCF2GDS(vcf.fn, "ccm4.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm4.gds")
ccm_pca <- snpgdsPCA(genofile, autosome.only= F)

ccm_pca$sample.id

sNames <- ccm_pca$sample.id
sNames <- sub(".DmagnaLRV01.filtered.sorted.nd.bam", "", sNames)

layout(matrix(c(1,2,3,4), 2, 2, byrow = T)) 
layout.show(n=4) 

plot(ccm_pca$eigenvect[,1], ccm_pca$eigenvect[,2], col=c(rep('red',10),rep('blue',10),rep('green',10)), pch=19, main = c("All"))

plot(ccm_pca$eigenvect[,1][c(1:10)], ccm_pca$eigenvect[,2][c(1:10)], col=c(rep('red',10)), pch=19, main = c("DG"))
#text(ccm_pca$eigenvect[,1][c(1:10)], ccm_pca$eigenvect[,2][c(1:10)], labels=sNames[c(1:10)])
plot(ccm_pca$eigenvect[,1][c(11:20)], ccm_pca$eigenvect[,2][c(11:20)], col=c(rep('blue',10)), pch=19, main = c("DM"))
#text(ccm_pca$eigenvect[,1][c(11:20)], ccm_pca$eigenvect[,2][c(11:20)], labels=sNames[c(11:20)])
plot(ccm_pca$eigenvect[,1][c(21:30)], ccm_pca$eigenvect[,2][c(21:30)], col=c(rep('green',10)), pch=19, main = c("DP"))
# text(ccm_pca$eigenvect[,1][c(21:30)], ccm_pca$eigenvect[,2][c(21:30)], labels=sNames[c(21:30)])


library(clusterProfiler)

rgsvp2 <- read.csv("C:\\Users\\ChangQing\\Desktop\\GO富集分析\\RGSVP2-WT.csv", header = T)
up_genes <- read.table("C:\\Users\\ChangQing\\Desktop\\GO富集分析\\Up_genes.txt",header = F)
up_genes <- as.character(up_genes$V1)
down_genes <- read.table("C:\\Users\\ChangQing\\Desktop\\GO富集分析\\Down_genes.txt",header = F)
down_genes <- as.character(down_genes$V1)

# AnnotationHub的用法，但这个不能做水稻的GO，因为LOC或RAP无法之间转换为ENTREZID
library(AnnotationHub)
library(AnnotationDbi)
hub <- AnnotationHub()
# save(hub, file = "Hub.Rdata")
query(hub, "oryza sativa")
rice <- hub[["AH72263"]]
#oryza_sativa <- rice
save(rice, file = "rice.Rdata")
length(keys(rice))[1]
columns(rice)
#library(AnnotationDbi)
#bitr(keys(rice)[1:10], 'ENTREZID', c("SYMBOL", "GO", "ONTOLOGY"), rice)
library(org.Osativa.eg.db) # 这个包不行

#up_genes_rap <- na.omit(up_genes_rap$RAP)
up_MF <- enrichGO(gene = up_genes,
                   OrgDb= org.Osativa.eg.db,
                   keyType = "GID",
                   pAdjustMethod = "none",
                   ont = "MF")
up_BP <- enrichGO(gene = up_genes,
                  OrgDb= org.Osativa.eg.db,
                  keyType = "GID",
                  pAdjustMethod = "none",
                  ont = "BP")
up_CC <- enrichGO(gene = up_genes,
                  OrgDb= org.Osativa.eg.db,
                  keyType = "GID",
                  ont = "CC")

dotplot(up_CC)
dotplot(up_BP)
dotplot(up_MF)

##################################################################
# 水稻MSU和RAP ID 转换在线网站：http://bioinf.mind.meiji.ac.jp/OryzaExpress/ID_converter.php
# 一种新方法
# 先用在线网站分析， http://systemsbiology.cau.edu.cn/agriGOv2/
# 再用ggplot绘图
up_genes_go <- read.table("C:\\Users\\ChangQing\\Desktop\\GO富集分析\\up_go_result.txt", sep = "\t", header = T)
down_genes_go <- read.table("C:\\Users\\ChangQing\\Desktop\\GO富集分析\\down_go_result.txt", sep = "\t", header = T)
library(ggplot2)
up_go <- ggplot(up_genes_go, aes(x=(Number.in.input.list/Number.in.BG.Ref),y=Description)) + geom_point(aes(size=Number.in.input.list,colour=Ontology))
up_go + labs(x="GeneRatio", y="", size="Count")
down_go <- ggplot(down_genes_go,aes(x=(Number.in.input.list/Number.in.BG.Ref),y=Description)) + geom_point(aes(size=Number.in.input.list,colour=Ontology))
down_go + labs(x="GeneRatio", y="", size="Count") 

#####################################
p <- ggplot(rgsv_mock_down[1:76,], aes(x=(Number.in.input.list/Number.in.BG.Ref),y=Description))+geom_point(aes(size=Number.in.input.list,colour=FDR))+scale_color_continuous(low="red",high="blue")
p+labs(x="Genetatio",y="BP",size="Count")
p1 <- ggplot(rgsv_down[1:74,], aes(x=(Number.in.input.list/Number.in.BG.Ref),y=Description))+geom_point(aes(size=Number.in.input.list,colour=FDR))+scale_color_continuous(low="red",high="blue")
p1+labs(x="Generatio",y="Biological Process",size="Count")
p2 <- ggplot(rgsv_up[1:71,], aes(x=(Number.in.input.list/Number.in.BG.Ref),y=Description))+geom_point(aes(size=Number.in.input.list,colour=FDR))+scale_color_continuous(low="red",high="blue")
p2+labs(x="Generatio",y="Biological Process",size="Count")

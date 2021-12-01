
rm(list = ls())
library(tidyverse)
library(lme4)
library(data.table)
library(lmerTest)
library(effectsize)
library(parallel)
library(pbmcapply)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

wd = "~/Documents/Phd/Thesis/ENIGMA/NEW_Analysis_MARCH2019/Virtual_histology/GeneCoExpression"
setwd(wd)
matrixdir="coexpr_matrix"
outdir = "grot_analysis"
if (!dir.exists(file.path("results",outdir))) { dir.create(file.path("results",outdir))}

options(bitmapType = "cairo")

squeue = function(user, job) {
  df.squeue = system(paste0("squeue --name=", job," -u ",user), intern = T) %>% 
  strsplit(., " +") %>% 
      simplify2array() %>% 
      t() %>% 
      as.data.frame()
  return(df.squeue)
  
}

```

## Preprocessing
```{r preproc}
if (!file.exists("data.wide.Rda")) {

# load data    
load("FullData/Expression_Data_ForYash3.rda")
df.consistentgenes = read_tsv("Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv")
df.consistentgenes$GeneSymbol = gsub("-", "_",df.consistentgenes$GeneSymbol)

# convert gene to character
Full_Data = 
  Full_Data %>% 
  mutate(Gene = as.character(Gene))

# remove genes without gene expression in all datasets
Full_Data.expected = 
  expand.grid(Gene = Full_Data$Gene %>% unique(),
              Database = Full_Data$Database %>% unique())

Full_Data.empirical = 
Full_Data %>% 
  group_by(Gene, Database) %>% 
  tally()

Full_data.incomplete.genes = 
  anti_join(Full_Data.expected, 
            Full_Data.empirical)

genes.remove = unique(Full_data.incomplete.genes$Gene) %>% as.character()

Full_Data = 
  Full_Data %>% 
  filter(!Gene %in% genes.remove)
Full_Data$Gene = gsub("-", "_",Full_Data$Gene)

# get unique genes and consistent genes
genes = Full_Data$Gene %>% unique() %>% as.character()
consistent.genes = genes[genes %in% df.consistentgenes$GeneSymbol]

# convert to wide format
Full_Data.wide = 
  Full_Data %>% 
  dplyr::select(-Expression) %>% 
  pivot_wider(names_from = Gene, 
              values_from = Scaled_Expression,
              values_fn = mean)

# save
save("Full_Data.wide", 
     "genes", 
     "consistent.genes",
     file = "data.wide.Rda")

# rm unnecessary files
rm(list = c("Full_Data.empirical", 
            "Full_Data.expected", 
            "Full_data.incomplete.genes", 
            "Full_Data"))
} else {
  load("data.wide.Rda") 
}

```

# coexpression
```{r coexpression}
df.out = data.frame(Gene = consistent.genes, 
                    stringsAsFactors = F)

proc_files = list.files(matrixdir) %>%  str_split(pattern = ".Rda", simplify = T) %>% .[,1]
proc_files = df.out$Gene[!df.out$Gene %in% proc_files]

if (length(proc_files) > 0 ) {
  pb <- txtProgressBar(min = 0, max = length(proc_files), style = 3)
  for (i in 1:length(proc_files)) {
    setTxtProgressBar(pb, i)
    system(paste("sbatch scripts/compute_coexpression.sh", 
                   proc_files[i],
                   wd, 
                   sep = " "))
      Sys.sleep(5)
  
      # check I am not overloading the nodes    
      df.squeue = squeue("p23-didacvp","compute_coexpression")
      while (length(df.squeue$V1) > 300) {
        Sys.sleep(120) 
        print("script running on sbatch")
        df.squeue = squeue("p23-didacvp","compute_coexpression")
      }
    }
}

```


# coexpression
```{r coexpression}

rank.coexpr <- function(cgene) {
  load(file.path("coexpr_matrix",
                 paste0(cgene, ".Rda")))
       
  tmp.coexpr = tmp %>% 
    filter(!is.na(T)) %>% 
    mutate(eta = t_to_eta2(T,df)[[1]], 
           eta = if_else(T > 0,eta,-eta),
           rank = rank(eta, ties.method = "average")/length(tmp$T)) %>% 
    filter(rank > .5) %>% 
    dplyr::select(Cgene, Tgene, rank)
    
  return(tmp.coexpr)
  }  

if (!file.exists("df.coexpression.Rda")) {

  df.coexpression = pbmclapply(df.out$Gene, function(x) {rank.coexpr(x)}, mc.cores = 10)
  df.coexpression = df.coexpression %>% 
    reduce(rbind) %>% 
    nest(cols = -Cgene)

  save(df.coexpression, file = "df.coexpression.Rda") 
  
} else {
  load("df.coexpression.Rda")  

}

```


# prepare phenotypic maps
# this will require adaptation.
```{r phenotype_map}
if (!file.exists(file.path("results",  outdir, "phenotype_data.Rda"))) { 
  df.deriv = read.table(file.path(wd, "der_thickness_specific_morethan100.csv"),
                 sep = ",",
                 header = T)
  
  df.group =
    data.frame(young = df.deriv %>% 
                 filter(X < 20) %>% 
                 dplyr::select(-X) %>% 
                 summarise_all(mean) %>% 
                 t(),
             middle = 
               df.deriv %>% 
               filter(X > 20 & X <60) %>% 
               dplyr::select(-X) %>% 
               summarise_all(mean) %>% 
               t(),
             old = df.deriv %>% 
               filter(X >60) %>% 
               dplyr::select(-X) %>% 
               summarise_all(mean) %>% 
               t(),
             regionames = 
               names(df.deriv)[-1] %>% 
               gsub("_thickavg", "",.) %>%  
               gsub("L_", "ctx-lh-", .)) %>% 
    pivot_longer(-regionames, 
                 names_to = "AgeGroup", 
                 values_to = "thickness") %>% 
    mutate(sign = if_else(AgeGroup == "old", "neg", "pos"))
  
  tmp = expand.grid(AgeGroup = unique(df.group$AgeGroup), 
              CellType = c("Microglia", "CA1.Pyramidal", "Astrocyte", "S1.Pyramidal", "Ependymal","Oligodendrocyte")) %>% 
    filter(!(AgeGroup == "young" & CellType %in% c("Ependymal", "Oligodendrocyte"))) %>%
        filter(!(AgeGroup == "middle" & CellType %in% c("Microglia", "CA1.Pyramidal"))) %>%
        filter(!(AgeGroup == "old" & CellType %in% c( "Astrocyte", "S1.Pyramidal", "Ependymal","Oligodendrocyte")))


  
  df.group = left_join(tmp, df.group) %>% 
    nest(data = c(regionames, thickness))
  df.group[df.group$AgeGroup=="middle" & (df.group$CellType=="Oligodendrocyte"|df.group$CellType=="S1.Pyramidal"),]$sign<-"neg"

  save(df.group, 
    file = file.path("results",
                     outdir,
                     "phenotype_data.Rda"))
} else {
  load(file.path("results",
                   outdir,
                   "phenotype_data.Rda"))
}
```


# Virtual Histology brain x gene data.
```{r VH data}
# 1 upload consistent genes and celltype
df.consistentgenes = read_tsv("Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv")
df.consistentgenes = df.consistentgenes %>% 
  mutate(GeneSymbol = gsub("-", "_",GeneSymbol)) %>% 
  filter(GeneSymbol %in% consistent.genes)


#1) et left hemisphere profiles of consistent genes.
tmp = read_tsv("AllenHBA_DK_ExpressionMatrix.tsv") 
tmp.GeneSymbol = tmp$X1 %>% gsub("-", "_",.)
tmp = tmp %>% dplyr::select(starts_with("ctx-lh")) %>% 
  t() %>% 
as.data.frame()
names(tmp) = tmp.GeneSymbol
df.profile.consistent.genes = tmp %>% 
  dplyr::select(consistent.genes)
```
```{r}
# 1 rank consistent genes depending on the profile. 
tmp.lj = data.frame(regionames = rownames(df.profile.consistent.genes))
df.spatial.correlation = list()
for (i in 1:length(df.group$AgeGroup)) {
  tmp = left_join(tmp.lj, df.group$data[[i]])
  cor.spatial = cor(tmp$thickness, df.profile.consistent.genes)
  
  tmp.df.spatial.correlation = 
    data.frame(
      AgeGroup = df.group$AgeGroup[i],
      CellType = df.group$CellType[i],
      sign = df.group$sign[i], 
      GeneSymbol = names(df.profile.consistent.genes),
      cor.spatial = cor.spatial %>% t())
  
  tmp.df.spatial.correlation = 
    tmp.df.spatial.correlation %>% mutate(rank.spatial = 
                                        if_else(sign == "pos", 
                                                rank(cor.spatial)/length(cor.spatial), 
                                                rank(-cor.spatial)/length(cor.spatial)))
 
  df.spatial.correlation[[i]] = tmp.df.spatial.correlation   
}
df.spatial.correlation = df.spatial.correlation %>% 
    purrr  ::reduce(rbind) %>% 
   nest(spatialcor = c(GeneSymbol,cor.spatial,rank.spatial))

df.group = left_join(df.group, df.spatial.correlation)
```

```{r overexpression}

thr.spatial = .95 # select genes that are 1) celltype specific and amongst the tope %5 associated with the phenotype of interest
thr.coexpression = .999 # select top .1% of genes coexpress with the genes of interest. 


get_genesCoexpressed = function(df.coexpression, df.consistentgenes,celltype, spatialcor, thr.spatial, thr.coexpression){
    genes.celltype = df.consistentgenes %>% 
    filter(CellType == celltype) %>% 
    .$GeneSymbol
  # top genes
  genes.top = 
    spatialcor %>% 
    filter(GeneSymbol %in% genes.celltype & rank.spatial > thr.spatial) %>% 
    .$GeneSymbol
  # top coexpressed genes
  genes.coexpressed = 
    df.coexpression %>% 
    filter(Cgene %in% genes.top) %>% 
    unnest() %>% 
    filter(rank > thr.coexpression) %>% 
    .$Tgene %>% 
    unique()
  return(genes.coexpressed)
}

# run over expression analysis 
co_exp_run <- function(genes.coexpressed, geneset, ngenes) {
  i <- length(intersect(genes.coexpressed, geneset))
  b <- length(geneset)
  d <- length(genes.coexpressed)
  out <- phyper(i-1, b, ngenes-b, d, lower.tail = F, log.p = F)
  return(out)
}

df.group$coexpressed.genes = NULL
for (i in 1:length(df.group$AgeGroup)) {
df.group$coexpressed.genes[[i]] = 
  get_genesCoexpressed(df.coexpression, 
                          df.consistentgenes,
                          celltype = df.group$CellType[i], 
                          spatialcor = df.group$spatialcor[[i]], 
                          thr.spatial, 
                          thr.coexpression)

}


save(df.group, 
     file = file.path("results",
                   outdir,
                   "ORAcoexpression.Rda"))

#clusterprofiler
#create a gene list for each age group
young=unique(unlist(df.group[df.group$AgeGroup=="middle",]$coexpressed.genes))
#getid

hs <- org.Hs.eg.db
my.symbols <- young
geneid<-AnnotationDbi::select(hs, 
               keys = my.symbols,
               columns = c("ENTREZID", "SYMBOL"),
               keytype = "SYMBOL")[,2]
v<-enrichGO(geneid, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01,pAdjustMethod="fdr",minGSSize=10,maxGSSize=500)
v2<-setReadable(v,org.Hs.eg.db)
dotplot(v)
cnetplot(v2)


####Different Thresholds

for (thr.spatial in seq(0.9,0.99,0.01)) {
for (thr.coexpression in seq(0.995,0.999,0.001)) {

for (i in 1:length(df.group$AgeGroup)) {
    genes.celltype = df.consistentgenes %>% 
    filter(CellType == df.group$CellType[i]) %>% 
    .$GeneSymbol
  # top genes
  genes.top = 
    df.group$spatialcor[[i]] %>% 
    filter(GeneSymbol %in% genes.celltype & rank.spatial > thr.spatial) %>% 
    .$GeneSymbol
  # top coexpressed genes
  if (length(genes.top)!=0) {genes.coexpressed = 
    df.coexpression %>% 
    filter(Cgene %in% genes.top) %>% 
    unnest() %>% 
    filter(rank > thr.coexpression) %>% 
    .$Tgene %>% 
    unique()
    if (length(genes.top)!=0) {df.group$coexpressed.genes[[i]] = 
genes.coexpressed}
else if  (length(genes.top)==0) {df.group$coexpressed.genes[[i]]=NA}
}
else {next}

}


temp<-df.group[!is.na(df.group$coexpressed.genes),]
for ( age in unique(temp$AgeGroup)) {
if (!is.null(unique(unlist(temp[temp$AgeGroup==age,]$coexpressed.genes))))

{young=unique(unlist(temp[temp$AgeGroup==age,]$coexpressed.genes))
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my.symbols <- young
geneid<-AnnotationDbi ::select(hs, 
               keys = my.symbols,
               columns = c("ENTREZID", "SYMBOL"),
               keytype = "SYMBOL")[,2]
v<-enrichGO(geneid, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01,pAdjustMethod="fdr",minGSSize=10,maxGSSize=500)
if (sum(v@result$p.adjust<0.01)!=0) {
v2<-setReadable(v,org.Hs.eg.db)
tiff(paste(thr.spatial,thr.coexpression, age, "BP_dotplot.tiff",sep=""),600,600)
print(dotplot(v,showCategory=20,title=paste(age,"-","thr.spatial=",thr.spatial,"thr.coexp=", thr.coexpression,sep=" ")))
dev.off()
}
else {next}

}
else {next}


}
}
}




# Each celltype and age group
temp<-df.group[!is.na(df.group$coexpressed.genes),]
for ( cell in unique(temp$CellType)) {
for ( age in unique(temp$AgeGroup)) {
if (!is.null(unique(unlist(temp[temp$AgeGroup==age & temp$CellType==cell,]$coexpressed.genes))))

{selected=unique(unlist(temp[temp$AgeGroup==age & temp$CellType==cell,]$coexpressed.genes))
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my.symbols <- selected
geneid<-AnnotationDbi ::select(hs, 
               keys = my.symbols,
               columns = c("ENTREZID", "SYMBOL"),
               keytype = "SYMBOL")[,2]
v<-enrichGO(geneid, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01,pAdjustMethod="fdr",minGSSize=10,maxGSSize=500)
if (sum(v@result$p.adjust<0.01)!=0) {
v2<-setReadable(v,org.Hs.eg.db)
tiff(paste(cell,thr.spatial,thr.coexpression, age, "BP_dotplot.tiff",sep=""),600,600)
print(dotplot(v,showCategory=20,title=paste(age,"-","thr.spatial=",thr.spatial,"thr.coexp", thr.coexpression, cell,sep=" ")))
dev.off()
}
else {next}

}
else {next}


}
}
}
}






#getid

```
```{r}
#sensitivity analysis- each celltype and age period with signficant
#cell="Astrocyte";lowage=4; highage=16; sign="pos"
#cell="CA1.Pyramidal";lowage=4; highage=16; sign="pos"
#cell="Microglia";lowage=4; highage=18; sign="pos"
#cell="S1.Pyramidal";lowage=4; highage=23; sign="neg"
#cell="Astrocyte";lowage=43; highage=61; sign="pos"
#cell="Ependymal";lowage=44; highage=61; sign="pos"
#cell="Oligodendrocyte";lowage=42; highage=46; sign="neg"
#cell="S1.Pyramidal";lowage=36; highage=64; sign="neg"
#cell="CA1.Pyramidal";lowage=68; highage=76; sign="neg"
#cell="Microglia";lowage=68; highage=76; sign="neg"



wd = "~/Documents/Phd/Thesis/ENIGMA/NEW_Analysis_MARCH2019/Virtual_histology/GeneCoExpression"
setwd(wd)
df.deriv = read.table(file.path(wd, "der_thickness_specific_morethan100.csv"),
                 sep = ",",
                 header = T)
  
  df.group =
    data.frame(young = df.deriv %>% 
                 filter(X < 20) %>% 
                 dplyr::select(-X) %>% 
                 summarise_all(mean) %>% 
                 t(),
             middle = 
               df.deriv %>% 
               filter(X > 20 & X <60) %>% 
               dplyr::select(-X) %>% 
               summarise_all(mean) %>% 
               t(),
             old = df.deriv %>% 
               filter(X >60) %>% 
               dplyr::select(-X) %>% 
               summarise_all(mean) %>% 
               t(),
             regionames = 
               names(df.deriv)[-1] %>% 
               gsub("_thickavg", "",.) %>%  
               gsub("L_", "ctx-lh-", .)) %>% 
    pivot_longer(-regionames, 
                 names_to = "AgeGroup", 
                 values_to = "thickness") %>% 
    mutate(sign = if_else(AgeGroup == "old", "neg", "pos"))
  
  tmp = expand.grid(AgeGroup = unique(df.group$AgeGroup), 
              CellType = cell)
  df.group = left_join(tmp, df.group) %>% 
    nest(data = c(regionames, thickness))


df.group<-df.group[1,]


df.group$data[[1]]$thickness<-rowMeans(t(df.deriv[df.deriv$X>=lowage & df.deriv$X<=highage,]))[-1]
df.group$sign=sign


df.consistentgenes = read_tsv("Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv")
df.consistentgenes = df.consistentgenes %>% 
  mutate(GeneSymbol = gsub("-", "_",GeneSymbol)) %>% 
  filter(GeneSymbol %in% consistent.genes)


#1) et left hemisphere profiles of consistent genes.
tmp = read_tsv("AllenHBA_DK_ExpressionMatrix.tsv") 
tmp.GeneSymbol = tmp$X1 %>% gsub("-", "_",.)
tmp = tmp %>% dplyr::select(starts_with("ctx-lh")) %>% 
  t() %>% 
as.data.frame()
names(tmp) = tmp.GeneSymbol
df.profile.consistent.genes = tmp %>% 
  dplyr::select(consistent.genes)
```
```{r}
# 1 rank consistent genes depending on the profile. 
tmp.lj = data.frame(regionames = rownames(df.profile.consistent.genes))
df.spatial.correlation = list()
for (i in 1:length(df.group$AgeGroup)) {
  tmp = left_join(tmp.lj, df.group$data[[i]])
  cor.spatial = cor(tmp$thickness, df.profile.consistent.genes)
  
  tmp.df.spatial.correlation = 
    data.frame(
      AgeGroup = df.group$AgeGroup[i],
      CellType = df.group$CellType[i],
      sign = df.group$sign[i], 
      GeneSymbol = names(df.profile.consistent.genes),
      cor.spatial = cor.spatial %>% t())
  
  tmp.df.spatial.correlation = 
    tmp.df.spatial.correlation %>% mutate(rank.spatial = 
                                        if_else(sign == "pos", 
                                                rank(cor.spatial)/length(cor.spatial), 
                                                rank(-cor.spatial)/length(cor.spatial)))
 
  df.spatial.correlation[[i]] = tmp.df.spatial.correlation   
}
df.spatial.correlation = df.spatial.correlation %>% 
    purrr  ::reduce(rbind) %>% 
   nest(spatialcor = c(GeneSymbol,cor.spatial,rank.spatial))

df.group = left_join(df.group, df.spatial.correlation)


thr.spatial = .95 # select genes that are 1) celltype specific and amongst the tope %5 associated with the phenotype of interest
thr.coexpression = .999 # select top .1% of genes coexpress with the genes of interest. 


get_genesCoexpressed = function(df.coexpression, df.consistentgenes,celltype, spatialcor, thr.spatial, thr.coexpression){
    genes.celltype = df.consistentgenes %>% 
    filter(CellType == celltype) %>% 
    .$GeneSymbol
  # top genes
  genes.top = 
    spatialcor %>% 
    filter(GeneSymbol %in% genes.celltype & rank.spatial > thr.spatial) %>% 
    .$GeneSymbol
  # top coexpressed genes
  genes.coexpressed = 
    df.coexpression %>% 
    filter(Cgene %in% genes.top) %>% 
    unnest() %>% 
    filter(rank > thr.coexpression) %>% 
    .$Tgene %>% 
    unique()
  return(genes.coexpressed)
}

# run over expression analysis 
co_exp_run <- function(genes.coexpressed, geneset, ngenes) {
  i <- length(intersect(genes.coexpressed, geneset))
  b <- length(geneset)
  d <- length(genes.coexpressed)
  out <- phyper(i-1, b, ngenes-b, d, lower.tail = F, log.p = F)
  return(out)
}

df.group$coexpressed.genes = NULL
for (i in 1:length(df.group$AgeGroup)) {
df.group$coexpressed.genes[[i]] = 
  get_genesCoexpressed(df.coexpression, 
                          df.consistentgenes,
                          celltype = df.group$CellType[i], 
                          spatialcor = df.group$spatialcor[[i]], 
                          thr.spatial, 
                          thr.coexpression)

}


#create a gene list for each age group
young=unique(unlist(df.group$coexpressed.genes))
#getid
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my.symbols <- young
geneid<-select(hs, 
               keys = my.symbols,
               columns = c("ENTREZID", "SYMBOL"),
               keytype = "SYMBOL")[,2]
v<-enrichGO(geneid, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01,pAdjustMethod="fdr",minGSSize=10,maxGSSize=500)

tiff(paste(cell,lowage,highage, "BP_dotplot.tiff",sep=""),600,600)
v2<-setReadable(v,org.Hs.eg.db)
print(dotplot(v,showCategory=20,title=paste(lowage,"-", highage,"years",cell,sep=" ")))
dev.off()

tiff(paste(cell,lowage,highage, "BP_cnetplot.tiff",sep=""),800,800)
print(cnetplot(v2,title=paste(lowage,"-", highage,"years",cell,sep=" ")))
dev.off()



v<-enrichGO(geneid, 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.01,pAdjustMethod="fdr",minGSSize=10,maxGSSize=500)

tiff(paste(cell,lowage,highage, "CC_dotplot.tiff",sep=""),600,600)
v2<-setReadable(v,org.Hs.eg.db)
print(dotplot(v,showCategory=20, title=paste(lowage,"-", highage,"years",cell,sep=" ")))
dev.off()

tiff(paste(cell,lowage,highage, "CC_cnetplot.tiff",sep=""),800,800)
print(cnetplot(v2,title=paste(lowage,"-", highage,"years",cell,sep=" ")))
dev.off()






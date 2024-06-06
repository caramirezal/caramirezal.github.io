## pseudotime analysis of brain organoids scRNA-Seq data
## https://www.nature.com/articles/s41586-022-05279-8

library(Biobase)

root.dir = '/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/'
data.dir = '/media/ag-cherrmann/cramirez'
out.dir = '/media/ag-cherrmann/cherrmann'

setwd(root.dir)

### read expression values

colon.exp.file = '/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/courses/data/colon_scenic_corrected_input.tsv'
colon.exp = read.table(colon.exp.file,header=TRUE,row.names=1)
colon.exp = t(colon.exp)

colon.meta.file = '/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/courses/data/colon_scenic_corrected_annotations.tsv'
colon.meta = read.table(colon.meta.file,header=TRUE,row.names=1)

table(colon.meta$orig.ident,colon.meta$Infection)

### read TF activities

colon.scenic.file = '/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/courses/data/colon_aucell.csv'
colon.scenic = read.csv(colon.scenic.file,row.names=1,check.names = FALSE)

### check that the order of cells is comparable

sum(rownames(colon.scenic) != rownames(colon.meta))  ## yes, they are!

colon.meta = data.frame(colon.meta,colon.scenic)
meta = AnnotatedDataFrame(colon.meta)


### create ExpressionSet

saveRDS(colon.exp, 
        "/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/courses/data/colon_exp.RDS",
        compress = TRUE
)

library(Biobase)
eset = ExpressionSet(assayData = colon.exp,
                     phenoData = meta)

#eset.sub = eset[,pData(eset)$orig.ident != 'Colon_Mock']  ## we remove the mock cells

#i.cell = which(pData(eset)$orig.ident == 'Colon_24h')
i.cell = which(pData(eset)$CellTypes == 'Inmature Enterocyte 2')

eset.sub = eset[,i.cell]  ## we remove the mock cells

### should we remove some weird genes?

m = apply(exprs(eset.sub),1,mean)
s = apply(exprs(eset.sub),1,sd)

plot(m,s,pch=20)
plot(density(m),lwd=3,col='red',type='l')





### Load the list of DEGs between infected and bystander cells

deg12.file = '/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/courses/data/degs_imm_infected_vs_bystander_12_colon.rds'
deg24.file = '/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/courses/data/degs_imm_infected_vs_bystander_24_colon.rds'

tmp12 = readRDS(deg12.file)
tmp24 = readRDS(deg24.file)

fdr = 0.05
deg12 = rownames(tmp12[tmp12$p_val_adj < fdr,])
deg24 = rownames(tmp24[tmp24$p_val_adj < fdr,])

deg = union(deg12,deg24)

### load the INF genes

ifn = read.table('/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/courses/data//REACTOME_INF_SIGNALLING.txt',as.is=TRUE)[,1]


### diffusion map

library(destiny)


# selecting the genes on which to compute the pseudo-time

i.top1 = which(s>2) # genes with a high variance
i.top2 = which(rownames(eset.sub) %in% deg)  # take the genes DEGs between infected and bystanders
i.top3 = which(rownames(eset.sub) %in% ifn) # take genes in the IFN pathway

length(i.top2)
length(i.top3)
length(union(i.top2,i.top3))

i.top = union(i.top1,i.top3)
i.top = i.top3

# compute the diffusion map
dm = DiffusionMap(eset.sub[i.top,],n_pcs = 50)


saveRDS(eset.sub, 
        "/Users/cbg-mbp-02/Documents/scripts/caramirezal.github.io/courses/data/sce_colon_immEnt.RDS", 
        compress = TRUE
)


plot(dm, 1:2,col_by='Infection')
plot(dm, 1:2,col_by='orig.ident')
plot(dm, 1:2,col_by='CellTypes')

plot(dm, 1:2,col_by='NFKB2...')
plot(dm, 1:2,col_by='STAT1...')
plot(dm, 1:2,col_by='IRF1...')
plot(dm, 1:2,col_by='JUN...')


plot(dm, 1:2,col_by='COVID19')

### compute diffusion pseudo-time
dpt = DPT(dm)

plot(dpt,col_by='DPT4')

### make a matrix with the DPT times

dpt.times = as.data.frame(dpt)
dim(dpt.times)

### correlate TF activities with DPTs
XXX = cor(colon.scenic[i.cell,],dpt.times[,1:50])

pheatmap(cor(t(XXX)),cex=.7)
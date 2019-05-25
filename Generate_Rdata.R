suppressMessages(require(data.table))
suppressMessages(require(Rtsne))
suppressMessages(require(Matrix))
suppressMessages(require(cellrangerRkit))
suppressMessages(require(feather))
suppressMessages(require(Seurat))
suppressMessages(require(matrixStats))
require(mixtools)
source("~/Analysis/201801_JohnVDJ/analysis/VDJ_function_pack.R")
devtools::load_all("/home/ahe/tools/Lightbulb")

#load("TSNE_loading.Rdata")
ptm <- proc.time()
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#parameter setting
genome <- "mm10"
cell_per_sample=2000
cellrangerOutDir="/home/ahe/Analysis/201702_John10X/data"
SampleIdPattern="^d[0-9]"
tag=dir(path =cellrangerOutDir,pattern=SampleIdPattern,full.names = F, recursive = F)
tagdir=paste0(cellrangerOutDir,"/",tag)[gtools::mixedorder(tag)]
nGeneThreshold=400
mitoUpperThreshold=0.2
mitoLowerThreshold=0.005
upper_ranking_threshold_pecentage=0.001
min_cellnum_per_cluster=21
##reference gene table
t2g=readRDS("/home/ahe/Analysis/genomeFiles/mm10_t2g_R35.Rds")

MTlist=t2g$ens_gene[t2g$chromosome_name=="MT"]
IGlist=t2g$ens_gene[grep("Igh|Igk|Igl.*-",t2g$ext_gene)]

#combining of all single library
set.seed(123)
gene_bc_matrix=list()
sample_num=rep(0,length(tagdir))
cell_passed_filtering=rep(0,length(sample_list))
cell_selected_for_downstream=rep(0,length(sample_list))
batch=c()
for(i in 1:length(tagdir)){
    #readin cells
    gene_bc_matrix[[i]]=load_cellranger_matrix_h5(tagdir[i], genome=genome,barcode_filtered =F)
    gene_bc_matrix[[i]]=exprs(gene_bc_matrix[[i]])
    gene_bc_matrix[[i]]=gene_bc_matrix[[i]][,colSums(gene_bc_matrix[[i]])>nGeneThreshold]
    #filter by nGene
    gene_bc_matrix[[i]]=t(gene_bc_matrix[[i]])
    temp=gene_bc_matrix[[i]][,colSums(gene_bc_matrix[[i]])>nrow(gene_bc_matrix[[i]])/100] # store cells with filtered gene num
    enoughGene=rowSums(temp>0)>=nGeneThreshold
    gene_bc_matrix[[i]]=gene_bc_matrix[[i]][enoughGene,] 
    MTsum=rowSums(gene_bc_matrix[[i]][,colnames(gene_bc_matrix[[i]]) %in% MTlist])
    othersum=rowSums(gene_bc_matrix[[i]][,!colnames(gene_bc_matrix[[i]]) %in% MTlist])
    toomuchmito=MTsum/rowSums(gene_bc_matrix[[i]])>mitoUpperThreshold
    toolittlemito=MTsum/othersum <= mitoLowerThreshold/(1-mitoLowerThreshold)
    gene_bc_matrix[[i]]=gene_bc_matrix[[i]][(!toomuchmito) & (!toolittlemito),]
    
    cell_passed_filtering[i]=nrow(gene_bc_matrix[[i]])
    if(nrow(gene_bc_matrix[[i]])>cell_per_sample){
        topick <- sample(1:nrow(gene_bc_matrix[[i]]),cell_per_sample)
        gene_bc_matrix[[i]]=gene_bc_matrix[[i]][topick,]
    }
    
    #some other statistics
    sample_num[i]=ncol(gene_bc_matrix[[i]]) #sample is now super cells
    batch=c(batch,rep(i,nrow(gene_bc_matrix[[i]])))
    rownames(gene_bc_matrix[[i]])=gsub("1$", i, rownames(gene_bc_matrix[[i]]))
    
    cell_selected_for_downstream[i]=nrow(gene_bc_matrix[[i]])
    #message    
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2)) 
    message(paste0(i,"/",length(tagdir)," loading done: ",t[1],"hr",t[2],"min",round(t[3]),"s"))
    message(paste0("droped by too few gene:",sum(!enoughGene),"| mito:",sum(toomuchmito),"|",sum(toolittlemito)))
    message(paste0("cell passed filtering: ",cell_passed_filtering[i]," | cell selected for downstream: ",cell_selected_for_downstream[i])," (",round(cell_selected_for_downstream[i]/cell_passed_filtering[i]*100,1),"%)")
}
message("")
message(paste0("in total: ",sum(cell_passed_filtering)," passed filtering & ",sum(cell_selected_for_downstream)," selected for downstream analysis"))

ExpressionMat=do.call(rbind,gene_bc_matrix)
batchidx=gsub("/home/ahe/Analysis/201702_John10X/data/","",tagdir)
rm(gene_bc_matrix)
ExpressionMat=ExpressionMat[,colSums(ExpressionMat)>2]
ExpressionBinaryMat=ExpressionMat>=1
ExpressionMat=ExpressionMat[,colSums(ExpressionBinaryMat)>nrow(ExpressionMat)/100]
rm(ExpressionBinaryMat)

#filter the outlier for each gene:
upper_ranking_threshold=ceiling(upper_ranking_threshold_pecentage*nrow(ExpressionMat))
ExpressionMat=apply(ExpressionMat,2,function(x){
    ordered_x=sort(x,decreasing=T)
    max_threshold=ordered_x[upper_ranking_threshold]
    x[x>max_threshold]=max_threshold
    return(x)
})

#change name
ExpressionMat=ExpressionMat[,which(colnames(ExpressionMat) %in% t2g$ens_gene)]
colnames(ExpressionMat)=t2g$ext_gene[match(colnames(ExpressionMat),t2g$ens_gene)]
colnames(ExpressionMat)[which(duplicated(colnames(ExpressionMat)))]=paste0(colnames(ExpressionMat)[which(duplicated(colnames(ExpressionMat)))],"_1")

#z score by column(gene)
ExpressionMat_TPM=t(apply(as.matrix(ExpressionMat),1,function(x){x/sum(x)*1000000}))

#generate seurat
Exp_Seurat <- CreateSeuratObject(raw.data = t(ExpressionMat))
rm(ExpressionMat)
Exp_Seurat@meta.data$batch <- batch
Exp_Seurat@data=t(log2(ExpressionMat_TPM+1))
Exp_Seurat@raw.data=as(Exp_Seurat@raw.data, "sparseMatrix") 
Exp_Seurat@data=as(Exp_Seurat@data, "sparseMatrix") 
Exp_Seurat@scale.data=t(scale(t(apply(as.matrix(Exp_Seurat@raw.data),2,function(x){x/sqrt(sum(x^2))*100000}))))
rm(ExpressionMat_TPM)
Exp_Seurat <- FindVariableGenes(Exp_Seurat, do.plot = F)
HVG <- head(rownames(Exp_Seurat@hvg.info), 5000)

ident_to_rm=c("placeholder")
while(length(ident_to_rm)>0){
    message(paste("round",i))
    #PCA
    ptm <- proc.time()
    Exp_Seurat <- RunPCA(object = Exp_Seurat, pc.genes = HVG[1:3000],pcs.compute = 100,do.print = F)
    #pca_allcell=prcomp(ExpressionNormed[,HVG[1:3000]])
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    message(paste("PCA done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

    Exp_Seurat <- FindClusters(object = Exp_Seurat, reduction.type = "pca", dims.use = 1:25,resolution = 1, print.output = F, save.SNN = TRUE, force.recalc =T)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    message(paste("Clustering done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

    #filter the singlets
    cluster_cell_num=table(Exp_Seurat@ident)
    cluster_names <- names(cluster_cell_num)
    ident_to_rm=cluster_names[cluster_cell_num<min_cellnum_per_cluster]
    if(length(ident_to_rm)>0){
        Exp_Seurat=filter_ident(Exp_Seurat,ident_to_rm)
        message(paste("Cleaning done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))
    }else{
        message("loop clean finished")
        break
    }
}

Exp_Seurat <- RunTSNE(object = Exp_Seurat, dims.use = 1:25, do.fast = TRUE)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("TSNE done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

Exp_Seurat <- umap_seurat(Exp_Seurat,pca_dim=25)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("umap done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

time_points=unique(gsub("_.*","",batchidx))
tissue_types=c("Spleen","Gut")

Exp_Seurat@meta.data$tissue="Spleen"
Exp_Seurat@meta.data$tissue[grep("TRM",batchidx[Exp_Seurat@meta.data$batch])]="Gut"
Exp_Seurat@meta.data$timepoint=batchidx[Exp_Seurat@meta.data$batch]
Exp_Seurat@meta.data$timepoint=gsub("_.*","",Exp_Seurat@meta.data$timepoint)


#MAGIC
forMagic="For_magic.feather"
afterMagic="Magic_out.feather"

temp_Exp=as.data.frame(t(as.matrix(Exp_Seurat@data)))
#magic
write_feather(temp_Exp, forMagic)
rm(temp_Exp)
system(paste0("python /home/ahe/tools/Lightbulb/src/Lightbulb_MagicWrapper.py --matx ",forMagic," --out ",afterMagic))
ExpressionMat_magic = data.matrix(read_feather(afterMagic))
colnames(ExpressionMat_magic)=gsub("MAGIC ", "", colnames(ExpressionMat_magic))
rownames(ExpressionMat_magic)=colnames(Exp_Seurat@data)

t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
message(paste("MAGIC done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))

lib2rm=setdiff(1:length(batchidx),grep("_2$",batchidx))
Exp_Seurat_rep=filter_metadata(Exp_Seurat,lib2rm,"batch")

lib2rm=grep("_2$",batchidx)
Exp_Seurat = filter_metadata(Exp_Seurat,lib2rm,"batch")

ExpressionMat_magic_rep=ExpressionMat_magic[match(colnames(Exp_Seurat_rep@data),rownames(ExpressionMat_magic)),]
ExpressionMat_magic=ExpressionMat_magic[match(colnames(Exp_Seurat@data),rownames(ExpressionMat_magic)),]

saveRDS(ExpressionMat_magic,file="ExpressionMat_magic.rds")
saveRDS(ExpressionMat_magic_rep,file="ExpressionMat_magic_rep.rds")

Exp_Seurat <- FindClusters(object = Exp_Seurat, reduction.type = "pca", dims.use = 1:30,resolution = 3, print.output = F, save.SNN = TRUE, force.recalc =T,reuse.SNN = F)

save(batchidx,Exp_Seurat,Exp_Seurat_rep,time_points,tissue_types,file="Exp_Seurat_all10_180917.Rdata")

Exp_Seurat@meta.data$GutORnot=as.numeric(Exp_Seurat@meta.data$tissue == "Gut" & Exp_Seurat@meta.data$timepoint %in% c("d4","d7","d10","d14","d21","d32","d60"))
Exp_Seurat@meta.data$SpleenORnot=as.numeric(Exp_Seurat@meta.data$tissue == "Spleen" & Exp_Seurat@meta.data$timepoint %in% c("d4","d7","d10","d14","d21","d32","d60"))

Exp_Seurat_Spleen = filter_metadata(Exp_Seurat,0,"SpleenORnot")
Exp_Seurat_Spleen <- FindVariableGenes(Exp_Seurat_Spleen, do.plot = F)
HVG <- head(rownames(Exp_Seurat_Spleen@hvg.info), 5000)
Exp_Seurat_Spleen <- RunPCA(object = Exp_Seurat_Spleen, pc.genes = HVG[1:3000],pcs.compute = 50,do.print = F)
Exp_Seurat_Spleen <- RunTSNE(object = Exp_Seurat_Spleen, dims.use = 1:25, do.fast = TRUE)
saveRDS(Exp_Seurat_Spleen,file="Exp_Seurat_Spleen.Rds")
rm(Exp_Seurat_Spleen)


Exp_Seurat_Gut = filter_metadata(Exp_Seurat,0,"GutORnot")
Exp_Seurat_Gut <- FindVariableGenes(Exp_Seurat_Gut, do.plot = F)
HVG <- head(rownames(Exp_Seurat_Gut@hvg.info), 5000)
Exp_Seurat_Gut <- RunPCA(object = Exp_Seurat_Gut, pc.genes = HVG[1:3000],pcs.compute = 50,do.print = F)
Exp_Seurat_Gut <- RunTSNE(object = Exp_Seurat_Gut, dims.use = 1:25, do.fast = TRUE)
saveRDS(Exp_Seurat_Gut,file="Exp_Seurat_Gut.Rds")
rm(Exp_Seurat_Gut)

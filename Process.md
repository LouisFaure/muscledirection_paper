Pagoda2\_analysis
================

  - [Load custom pipeline functions and
    libraries](#load-custom-pipeline-functions-and-libraries)
  - [Load data](#load-data)
  - [Filter and process data using
    pagoda2](#filter-and-process-data-using-pagoda2)
  - [Pagoda2 pathway overdispersion
    analysis](#pagoda2-pathway-overdispersion-analysis)
  - [CytoTRACE analysis](#cytotrace-analysis)

# Load custom pipeline functions and libraries

``` r
library(pagoda2)
library(Seurat)
library(ggthemes)
doUMAP <- function(PCA,n_neighbors,min_dist,max_dim=2,seed.use=42){
  require(reticulate)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
    py_set_seed(seed = seed.use)
  }
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(n_neighbors = as.integer(x = n_neighbors), 
                           n_components = as.integer(x = max_dim), metric = "correlation", 
                           min_dist = min_dist)
  
  umap_output <- umap$fit_transform(as.matrix(x = PCA))
  rownames(umap_output)=rownames(PCA)
  colnames(umap_output)=paste0("UMAP",1:max_dim)
  
  return(umap_output)
}

p2wrapper <- function(counts,n_neighbors=30,min_dist=0.3,k=100,npcs=200,selpc=T,...) {
  rownames(counts) <- make.unique(rownames(counts))
  p2 <- Pagoda2$new(counts,log.scale=FALSE,n.cores=parallel::detectCores()/2,...)
  p2$adjustVariance(plot=T,gam.k=10)
  p2$calculatePcaReduction(nPcs=npcs,n.odgenes=NULL,maxit=1000)
  if (selpc){
    x <- cbind(1:npcs, p2$misc$PCA$d)
    line <- x[c(1, nrow(x)),]
    proj <- princurve::project_to_curve(x, line)
    optpoint <- which.max(proj$dist_ind)-1
    dev.new(width=5, height=4)
    par(mfrow=c(1,1))
    plot(x,xlab="PC", ylab="Variance explained")
    abline(v=optpoint,lty=2,col=2)
    cat(paste0(optpoint," PCs retained\n"))
    npcs=optpoint
    p2$calculatePcaReduction(use.odgenes = T, name='PCA', 
                             nPcs=optpoint, maxit=1000)
  }
  
  p2$makeKnnGraph(k=k,type='PCA',center=T,distance='cosine');
  p2$getKnnClusters(method=conos::leiden.community,type='PCA',name = "leiden",resolution=.5)
  
  # Produce UMAP embedding
  cat("Computing UMAP... ")
  p2$embeddings$PCA$UMAP=doUMAP(PCA = p2$reductions$PCA[,1:npcs],n_neighbors = n_neighbors,min_dist = min_dist)
  cat("done\n")
  invisible(p2)
}
```

# Load data

``` r
pka100=Seurat::Read10X_h5("100pka_cellbender_filtered.h5")
pka05=Seurat::Read10X_h5("05pka_cellbender_filtered.h5")


colnames(pka100)=paste0("hard_",colnames(pka100))
colnames(pka05)=paste0("soft_",colnames(pka05))
```

# Filter and process data using pagoda2

``` r
cd <- lapply(list(pka100,pka05),function(y) {
  y=gene.vs.molecule.cell.filter(y,max.cell.size = 5e5,min.cell.size = 1e3)
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = y), value = TRUE)
  percent.mito <- Matrix::colSums(y[mito.genes, ])/Matrix::colSums(y)
  y = y[,percent.mito<.1]
  return(y)})

counts <- do.call(cbind,cd)

p2l=lapply(cd,function(x) p2wrapper(as(x,"dgCMatrix")))

batch=factor(sapply(strsplit(colnames(counts),"_"),"[[",1))
names(batch)=colnames(counts)

p2=p2wrapper(counts,batch=batch)


#p2w=p2webwrapper(p2)
adata <- SingleCellExperiment(list(counts=t(p2$misc$rawCounts)),
    colData=DataFrame(leiden=as.character(p2$clusters$PCA$leiden),
                      condition=sapply(strsplit(rownames(p2$counts),"_"),"[[",1)),
    metadata=list(study="combined")
)


reducedDims(adata) <- list(pca=p2$reductions$PCA,
                         umap=p2$embeddings$PCA$UMAP)

batch=adata$condition;names(batch)=rownames(p2$counts);batch=factor(batch)

save(adata,file="_Output/adata_combined.RData")
```

# Pagoda2 pathway overdispersion analysis

``` r
p2$n.cores=10
go.env <- p2.generate.mouse.go(p2)
p2$testPathwayOverdispersion(setenv = go.env,
                                     recalculate.pca=F,
                                     correlation.distance.threshold = 0.95,verbose = T)

myGeneNames <- colnames(p2$counts)
goSets <- p2.generate.mouse.go.web(myGeneNames)
deSets <- get.de.geneset(p2, groups = p2$clusters$PCA$leiden, prefix = 'de_')
geneSets <- c(goSets, deSets)

palette=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9")

additionalMetadata=list()
additionalMetadata$leiden <- p2.metadata.from.factor(p2$clusters$PCA$leiden, displayname = 'Leiden', s = 0.7, v = 0.8,start = 0, end = 0.5,
                                                     pal = palette)

additionalMetadata$condition <- p2.metadata.from.factor(batch, 
                                                  displayname = 'Condition', 
                                                  s = 0.7, 
                                                  v = 0.8,
                                                  start = 0, 
                                                  end = 0.5)


p2w <- make.p2.app(
  p2,
  dendrogramCellGroups = p2$clusters$PCA$leiden,
  additionalMetadata = additionalMetadata,
  geneSets = geneSets,
  show.clusters = FALSE # Hide the clusters that were used for the dendrogram from the metadata
)


p2w$serializeToStaticFast("_Output/p2w_batch_corrected.bin")
save(p2w,file="_Output/p2w_batch_corrected.RData")
```

# CytoTRACE analysis

``` r
p2<-p2l[[1]]
results = CytoTRACE(t(as.matrix(p2$misc$rawCounts)))
plotCytoTRACE(results,emb=p2$embeddings$PCA$UMAP,outputDir = "figures/100pka")


adata <- SingleCellExperiment(list(counts=t(p2$misc$rawCounts)),
    colData=DataFrame(leiden=as.character(p2$clusters$PCA$leiden_a[rownames(p2$counts)])),
    metadata=list(study="100pka")
)

reducedDims(adata) <- list(pca=p2$reductions$PCA,
                         umap=p2$embeddings$PCA$UMAP)

save(adata,file="_Output/adata_100pka.RData")

p2<-p2l[[2]]
results = CytoTRACE(t(as.matrix(p2$misc$rawCounts)))
plotCytoTRACE(results,emb=p2$embeddings$PCA$UMAP,outputDir = "figures/05pka")

adata <- SingleCellExperiment(list(logcounts=t(p2$counts)),
    colData=DataFrame(leiden=as.character(p2$clusters$PCA$leiden_a[rownames(p2$counts)])),
    metadata=list(study="05pka")
)

reducedDims(adata) <- list(pca=p2$reductions$PCA,
                         umap=p2$embeddings$PCA$UMAP)

save(adata,file="_Output/adata_05pka.RData")
```

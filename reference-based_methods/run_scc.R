suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(singleCellNet))
library(scClassify)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ramify))

read_h5<-function(h5file){
    library(reticulate)
    handle <- import("anndata", convert = FALSE)
    ad <- handle$read_h5ad(h5file)

    obs <- py_to_r(ad$obs)
    if ("n_counts" %in% colnames(x = obs)) {
    colnames(x = obs) <- gsub(
          pattern = "n_counts",
          replacement = "nUMI",
          x = colnames(x = obs)
        )
    }

    genenames <- rownames(py_to_r(ad$var))
    cellnames <- rownames(obs)

    ncells = length(cellnames)
    ngenes = length(genenames)
    # Read X
    if("indptr" %in% names(ad$X)){
        ad.X <- sparseMatrix(
            i = as.numeric(x = ad$X$indices),
            p = as.numeric(x = ad$X$indptr),
            x = as.numeric(x = ad$X$data),
            index1 = FALSE,
            dims=c(ngenes,
                   ncells)
        )
    }else{
        ad.X <- t(py_to_r(ad$X))
    }
    rownames(ad.X) <- genenames
    colnames(ad.X) <- cellnames

    # Read raw
    ad.raw.X = NA
    if("raw" %in% names(ad)){
        if("X" %in% names(ad$raw)){
            genenames <- rownames(py_to_r(ad$raw$var))
            ngenes = length(genenames)
            if("indptr" %in% names(ad$raw$X)){
                ad.raw.X <- sparseMatrix(
                    i = as.numeric(x = ad$raw$X$indices),
                    p = as.numeric(x = ad$raw$X$indptr),
                    x = as.numeric(x = ad$raw$X$data),
                    index1 = FALSE,
                    dims=c(ngenes,
                           ncells)
                )
            }else{
                ad.raw.X <- t(py_to_r(ad$raw$X))
            }
            rownames(ad.raw.X) <- genenames
            colnames(ad.raw.X) <- cellnames
        }else{
            genenames <- rownames(py_to_r(ad$var))
            ngenes = length(genenames)
            if("indptr" %in% names(ad$X)){
                ad.raw.X <- sparseMatrix(
                    i = as.numeric(x = ad$X$indices),
                    p = as.numeric(x = ad$X$indptr),
                    x = as.numeric(x = ad$X$data),
                    index1 = FALSE,
                    dims=c(ngenes,
                           ncells)
                )
            }else{
                ad.raw.X <- t(py_to_r(ad$X))
            }
            rownames(ad.raw.X) <- genenames
            colnames(ad.raw.X) <- cellnames
        }
    }
    return(list(X=ad.X, rawX=ad.raw.X, obs=obs))
}

run_scClassify<-function(obsRef,expRef,obsPred,expPred){
    expRef <- as(expRef, "dgCMatrix")
    expPred <- as(expPred, "dgCMatrix")
    commonGenes = intersect(rownames(expRef), rownames(expPred))
    print(sprintf("number of common genes: %d",length(commonGenes)))
    expRef = expRef[commonGenes,]
    obsRef<-data.frame(cell=obsRef, x = obsRef)
    rownames(obsRef)<-colnames(expRef)
    obsRef<-droplevels(obsRef)
    obsList <- splitCommon(sampTab=obsRef, ncells=100, dLevel="cell")
    obsTrain = obsList[[1]]
    expTrain = expRef[,rownames(obsTrain)]
    t0<-system.time(scClassify_res <- scClassify(exprsMat_train = expTrain,
                      cellTypes_train = obsTrain$cell,
                      exprsMat_test = list(Pred = expPred),
                      cellTypes_test = list(Pred = obsPred),
                      tree = "HOPACH",
                      algorithm = "WKNN",
                      selectFeatures = c("limma"),
                      similarity = c("pearson"),
                      returnList = FALSE,
                      verbose = FALSE))
    print(t0)
    return(scClassify_res$testRes$Pred$pearson_WKNN_limma$predRes)
}

run_singleCellNet<-function(obsRef,expRef,expPred){
    expRef <- as(expRef, "dgCMatrix")
    expPred <- as(expPred, "dgCMatrix")
    commonGenes = intersect(rownames(expRef), rownames(expPred))
    print(sprintf("number of common genes: %d",length(commonGenes)))
    expRef = expRef[commonGenes,]
    obsRef<-data.frame(cell=obsRef, x = obsRef)
    rownames(obsRef)<-colnames(expRef)
    obsRef<-droplevels(obsRef)
    obsList = splitCommon(sampTab=obsRef, ncells=100, dLevel="cell")
    obsTrain = obsList[[1]]
    expTrain = expRef[,rownames(obsTrain)]

    t0<-system.time(class_info<-scn_train(stTrain = obsTrain, expTrain = expTrain, nTopGenes = 10, 
                                      nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "cell"))
    print(t0)
    #validate data
    # obsTestList = splitCommon(sampTab=obsList[[2]], ncells=100, dLevel="cell") #normalize validation data so that the assessment is as fair as possible
    # obsTest = obsTestList[[1]]
    # expTest = expRef[commonGenes,rownames(obsTest)]
    #predict
    classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expPred, nrand = 50)

    return(classRes_val_all)
}

suppressPackageStartupMessages(library(SingleR))
library(SingleCellExperiment)
run_SingleR<-function(obsRef,expRef,obsPred,expPred){
    expRef <- as(expRef, "dgCMatrix")
    expPred <- as(expPred, "dgCMatrix")
    commonGenes = intersect(rownames(expRef), rownames(expPred))
    print(sprintf("number of common genes: %d",length(commonGenes)))
    expRef = expRef[commonGenes,]
    obsRef<-data.frame(cell=obsRef, x = obsRef)
    rownames(obsRef)<-colnames(expRef)
    obsRef<-droplevels(obsRef)
    obsList <- splitCommon(sampTab=obsRef, ncells=100, dLevel="cell")
    obsTrain = obsList[[1]]
    expTrain = expRef[,rownames(obsTrain)]
    expPred = expPred[commonGenes,]

    sceTrain <- SingleCellExperiment(
            assays = list(
                counts = expm1(expTrain),
                logcounts = expTrain
            ), 
            colData = obsTrain
        )

    scePred <- SingleCellExperiment(
            assays = list(
                counts = expm1(expPred),
                logcounts = expPred
            ), 
            colData = obsPred
        )

    t0<-system.time(prediction <- SingleR(test = scePred, ref = sceTrain, labels = obsTrain$cell))
    print(t0)
    return(prediction)
}

suppressPackageStartupMessages(library(CHETAH))
run_CHETAH<-function(obsRef,expRef,expPred){
    expRef <- as(expRef, "dgCMatrix")
    expPred <- as(expPred, "dgCMatrix")
    commonGenes = intersect(rownames(expRef), rownames(expPred))
    print(sprintf("number of common genes: %d",length(commonGenes)))
    expRef = expRef[commonGenes,]
    obsRef<-data.frame(cell=obsRef, x = obsRef)
    rownames(obsRef)<-colnames(expRef)
    obsRef<-droplevels(obsRef)
    obsList = splitCommon(sampTab=obsRef, ncells=100, dLevel="cell")
    obsTrain = obsList[[1]]
    expTrain = expRef[,rownames(obsTrain)]
    expPred = expPred[commonGenes,]

    sceTrain <- SingleCellExperiment(
            assays = list(
                counts = expm1(expTrain),
                logcounts = expTrain
            ), 
            colData = DataFrame(celltypes = obsTrain$cell)
        )

    scePred <- SingleCellExperiment(
            assays = list(
                counts = expm1(expPred),
                logcounts = expPred
            ), 
            colData = obsPred
        )
    t0<-system.time(scePred <- CHETAHclassifier(input = scePred, ref_cells = sceTrain))
    print(t0)
    return(scePred$celltype_CHETAH)
}
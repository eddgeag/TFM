


library(reshape2)
library(ggplot2)
library(dplyr)
library(vegan)
library(limma)
library(ggplot2)
library(ggpubr)
ruta_directorio <-
  "../../datos/resultados_analisis_datos_procesados"
ruta_archivo_datos <- "datos_integromics_procesados.rds"
ruta_actual_directorio <- "./12-2-22/scripts_analysis_reunion"

datos <- readRDS(file.path(ruta_directorio, ruta_archivo_datos))

datos.genero <- datos$datos.genero
datos.genero.abs <- datos$datos.genero.abs
datos.filo <- datos$datos.filo
datos.filo.abs <- datos$datos.filo.abs
grupo <- datos$grupo.common
obesidad <- datos$obesidad.common

analisis_microbiota <-
  function(microbiota,
           grupo,
           obesidad,
           microbiota.abs) {
    mc.instances <- numMCInstances(microbiota)
    mc.all <- getMonteCarloInstances(microbiota)
    res <- vector("list", length = mc.instances)
    res2 <- vector("list", length = mc.instances)
    res3 <- vector("list", length = mc.instances)
    res4 <- vector("list", length = mc.instances)
    res5 <- vector("list", length = mc.instances)
    
    

    matrix.obesity <- matrix(NA,ncol = 5,nrow=mc.instances)
    matrix.grupo <- matrix(NA,ncol = 5,nrow=mc.instances)
    matrix.interaccion <- matrix(NA,ncol = 5,nrow=mc.instances)

    
    
    for (mc.i in 1:mc.instances) {
      print(mc.i)
      t.input <- t(sapply(mc.all, function(y) {
        y[, mc.i]
      }))
      
      
      t.input.v <- (t(t.input - min(t.input)))
      t.input.none <- t(voom(t.input.v, normalize.method = "none")$E)
      t.input.scale <- t(voom(t.input.v, normalize.method = "scale")$E)
      t.input.quantile <- t(voom(t.input.v, normalize.method = "quantile")$E)
      t.input.cyclicloess <- t(voom(t.input.v, normalize.method = "cyclicloess")$E)
      
      
      
      
      mano <-
        adonis2(t.input ~ grupo * obesidad, method = "euclidean")
      mano2 <-
        adonis2(t.input.none ~ grupo * obesidad, method = "euclidean")
      mano3 <-
        adonis2(t.input.scale ~ grupo * obesidad, method = "euclidean")
      mano4 <-
        adonis2(t.input.quantile ~ grupo * obesidad, method = "euclidean")
      mano5 <-
        adonis2(t.input.cyclicloess ~ grupo * obesidad, method = "euclidean")
      
      matrix.obesity[mc.i,1] <-mano$`Pr(>F)`[2]
      matrix.obesity[mc.i,2] <-mano2$`Pr(>F)`[2]
      matrix.obesity[mc.i,3] <-mano3$`Pr(>F)`[2]
      matrix.obesity[mc.i,4] <-mano4$`Pr(>F)`[2]
      matrix.obesity[mc.i,5] <-mano5$`Pr(>F)`[2]
      
      matrix.grupo[mc.i,1] <-mano$`Pr(>F)`[1]
      matrix.grupo[mc.i,2] <-mano2$`Pr(>F)`[1]
      matrix.grupo[mc.i,3] <-mano3$`Pr(>F)`[1]
      matrix.grupo[mc.i,4] <-mano4$`Pr(>F)`[1]
      matrix.grupo[mc.i,5] <-mano5$`Pr(>F)`[1]
      
      matrix.interaccion[mc.i,1] <-mano$`Pr(>F)`[3]
      matrix.interaccion[mc.i,2] <-mano2$`Pr(>F)`[3]
      matrix.interaccion[mc.i,3] <-mano3$`Pr(>F)`[3]
      matrix.interaccion[mc.i,4] <-mano4$`Pr(>F)`[3]
      matrix.interaccion[mc.i,5] <-mano5$`Pr(>F)`[3]
      
   
      
      
      res[[mc.i]] <- log2(t.input-min(t.input)+2)
      res2[[mc.i]] <- t.input.none
      res3[[mc.i]] <- t.input.scale
      res4[[mc.i]] <- t.input.quantile
      res5[[mc.i]] <- t.input.cyclicloess
      
      
      
      
    }
    
    
    medias <- microbiota.abs * 0
    
    aux <- Reduce("rbind", res)
    for (j in 1:ncol(medias)) {
      medias[, j] <-
        rowMeans(matrix(aux[, j], nrow = nrow(medias), ncol = mc.instances))
      
      
    }
    colnames(medias) <- colnames(microbiota.abs)
    rownames(medias) <- rownames(microbiota.abs)
    
    medias2 <- microbiota.abs * 0
    
    aux <- Reduce("rbind", res2)
    for (j in 1:ncol(medias)) {
      medias2[, j] <-
        rowMeans(matrix(aux[, j], nrow = nrow(medias), ncol = mc.instances))
      
      
    }
    colnames(medias2) <- colnames(microbiota.abs)
    rownames(medias2) <- rownames(microbiota.abs)
    
    medias3 <- microbiota.abs * 0
    
    aux <- Reduce("rbind", res3)
    for (j in 1:ncol(medias)) {
      medias3[, j] <-
        rowMeans(matrix(aux[, j], nrow = nrow(medias), ncol = mc.instances))
      
      
    }
    colnames(medias3) <- colnames(microbiota.abs)
    rownames(medias3) <- rownames(microbiota.abs)
    
    medias4 <- microbiota.abs * 0
    
    aux <- Reduce("rbind", res4)
    for (j in 1:ncol(medias)) {
      medias4[, j] <-
        rowMeans(matrix(aux[, j], nrow = nrow(medias), ncol = mc.instances))
      
      
    }
    colnames(medias4) <- colnames(microbiota.abs)
    rownames(medias4) <- rownames(microbiota.abs)
    
    medias5 <- microbiota.abs * 0
    
    aux <- Reduce("rbind", res5)
    for (j in 1:ncol(medias)) {
      medias5[, j] <-
        rowMeans(matrix(aux[, j], nrow = nrow(medias), ncol = mc.instances))
      
      
    }
    colnames(medias5) <- colnames(microbiota.abs)
    rownames(medias5) <- rownames(microbiota.abs)

    adonis.obesidad <- colMeans(matrix.obesity)
    adonis.grupo <- colMeans(matrix.grupo)
    adonis.interaccion <- colMeans(matrix.interaccion)
    
    multivariante <- list(
      obesidad = adonis.obesidad,
      grupo = adonis.grupo,
      interaccion = adonis.interaccion,
      novoom = medias,
      voomnone=medias2,
      voomscale=medias3,
      voomquantile= medias4,
      voomcyclic = medias5
    )
    return(multivariante)
  }


res <-
  analisis_microbiota(
    microbiota = datos.genero,
    grupo,
    obesidad,
    microbiota.abs = datos.genero.abs
  )



ploteo <- function(csa){

  pcx <- prcomp(t(csa))
  X <- bind_cols(pcx$rotation[,1:2],obesidad=obesidad,grupo=grupo)
  csa$grupo <- grupo
  csa$obesidad <- obesidad
  csa.melt <- melt(csa)
  p1<-ggplot(csa.melt,aes(x=value,color=grupo))+geom_density()+facet_wrap(~obesidad)
  p2 <-ggplot(csa.melt,aes(y=value,color=grupo))+geom_boxplot()+facet_wrap(~obesidad)
  p3 <- ggplot(X,aes(PC1,PC2,color=obesidad,shape=grupo))+geom_point()
  ggarrange(p1,p2,p3)
}


ploteo((res$novoom))
ploteo(log2(res$voomnone+7))
ploteo(res$voomscale)
ploteo(res$voomquantile)
ploteo(res$voomcyclic)


















library(ALDEx2)
library(vegan)
library(reshape2)
library(ggplot2)
library(limma)
library(ggpubr)
construct_pre_Results <- function(X, correct, grupo, obesidad) {
  resultado <-
    inferencia_univariante(X,
                           grupo = grupo,
                           obesidad = obesidad,
                           progress = F)
  p.values <- resultado$p.values
  diferencias <- resultado$diferencias
  tipos_de_pruebas <- resultado$tipo_prueba
  if (!correct) {
    
  } else{
    p.values[, 1:4] <-
      apply(p.values[, 1:5], 2, function(x)
        p.adjust(x, "BH"))
    p.values[, 5:ncol(p.values)] <-
      apply(p.values[, 5:ncol(p.values)], 1, function(x)
        p.adjust(x, "BH"))
    
  }
  modelo <- bind_cols(tipos_de_pruebas, diferencias, p.values)
  
  return(modelo)
  
}

ruta_directorio <-
  "../../datos/resultados_analisis_datos_procesados"
ruta_archivo_datos <- "datos_integromics_procesados.rds"
ruta_actual_directorio <- "./12-2-22/scripts_analysis_reunion"
funciones_inferencia <- "infererencia_univariante_funcion.R"
source(file.path(ruta_actual_directorio, funciones_inferencia))
datos <- readRDS(file.path(ruta_directorio, ruta_archivo_datos))
grupo <- datos$grupo.common
obesidad <- datos$obesidad.common

analisis_microbiota <-
  function(microbiota,
           grupo,
           obesidad,
           vom,
           method,
           correccion) {
    
    microbiota<-datos$datos.genero
    vom<-F
    correccion<-F

    mc.instances <- numMCInstances(microbiota)
    mc.all <- getMonteCarloInstances(microbiota)
    res <- vector("list", length = length(names(microbiota@reads)))
    # res <- vector("list", length = 2)
    
    vector.obesity <-
      vector("numeric", length = length(names(microbiota@reads)))
    vector.grupo <-
      vector("numeric", length = length(names(microbiota@reads)))
    vector.interaccion <-
      vector("numeric", length = length(names(microbiota@reads)))
    
    
    for (mc.i in 1:mc.instances) {
      print(mc.i)
      t.input <- (sapply(mc.all, function(y) {
        y[, mc.i]
      }))
      
      if (vom) {
        t.input <- t(t.input - min(t.input))
        t.input <- t(voom(t.input, normalize.method = method)$E)
        
        
      }
      
      mano <-
        adonis2(t.input ~ grupo * obesidad, method = "euclidean")
      vector.grupo[mc.i] <- mano$`Pr(>F)`[1]
      vector.obesity[mc.i] <- mano$`Pr(>F)`[2]
      vector.interaccion[mc.i] <- mano$`Pr(>F)`[3]
      
      res[[mc.i]] <-
        construct_pre_Results(t.input,
                              correct = correccion,
                              grupo = grupo,
                              obesidad = obesidad)
      
      
      
    }
    
    
    diferencias <- lapply(res, function(x)
      x[15:27])
    pvalores <- lapply(res, function(x)
      x[28:44])
    
    pvaloresaux <- t(apply(Reduce("cbind", pvalores),1,as.numeric))
    diferenciasaux <-
      t(apply(Reduce("cbind", diferencias),1,as.numeric))
    
 
    matriz_pval <-
      matrix(NA,
             nrow = ncol(microbiota@reads),
             ncol = ncol(pvalores[[1]]))
    matriz_diferencias <-
      matrix(NA,
             nrow = ncol(microbiota@reads),
             ncol = ncol(diferencias[[1]]))
    
    nombres <- colnames(pvalores[[1]])
    nombresdif <- colnames(diferencias[[1]])
    
    for (col in 1:ncol(pvalores[[1]])) {
      nombre_i <- nombres[col]
      matriz_pval[, col] <-
        rowMeans(pvaloresaux[, grep(nombre_i, colnames(pvaloresaux))])
      
    }
    rownames(matriz_pval) <- names(microbiota@reads)
    colnames(matriz_pval) <- colnames(pvalores[[1]])
    
    for (col in 1:ncol(diferencias[[1]])) {
      nombresdif_i <- nombresdif[col]
      matriz_diferencias[, col] <-
        rowMeans(diferenciasaux[, grep(nombresdif_i, x = colnames(diferenciasaux))])
      
    }
    
    rownames(matriz_diferencias) <- names(microbiota@reads)
    colnames(matriz_diferencias) <- colnames(diferencias[[1]])
    adonis.obesidad <- mean(vector.obesity)
    adonis.grupo <- mean(vector.grupo)
    adonis.interaccion <- mean(vector.interaccion)
    
    retorno <-
      list(
        pvalues = matriz_pval,
        matriz_diferencias = matriz_diferencias,
        adonis.interaccion = adonis.interaccion,
        adonis.grupo = adonis.grupo,
        adonis.interaccion = adonis.interaccion
      )
  }


microbiota <- datos$datos.filo

res1_no_voom_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = F,
    method = "none",
    correccion = F
  )
print("######### PRIMERA ########")
res1_no_voom_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = F,
    method = "none",
    correccion = T
  )
print("######### SEGUNDA ########")

res1_voom_none_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "none",
    correccion = F
  )
print("######### TERCERA ########")

res1_voom_none_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "none",
    correccion = T
  )

print("######### CUARTA ########")

res1_voom_scale_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "scale",
    correccion = F
  )
print("######### QUINTA ########")

res1_voom_scale_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "scale",
    correccion = T
  )
print("######### SEXTA ########")

res1_voom_quantile_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "quantile",
    correccion = F
  )

res1_voom_quantile_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "quantile",
    correccion = T
  )
print("######### SEPTIMA ########")

res1_voom_quantile_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "cyclicloess",
    correccion = F
  )
print("######### OCTAVA ########")

res1_voom_none_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "cyclicloess",
    correccion = T
  )

microbiota <- datos$datos.genero

res1_no_voom_none_genero <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = F,
    method = "none",
    correccion = F
  )
print("######### PRIMERA ########")
res1_no_voom_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = F,
    method = "none",
    correccion = T
  )
print("######### SEGUNDA ########")

res1_voom_none_none_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "none",
    correccion = F
  )
print("######### TERCERA ########")

res1_voom_none_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "none",
    correccion = T
  )

print("######### CUARTA ########")

res1_voom_scale_none_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "scale",
    correccion = F
  )
print("######### QUINTA ########")

res1_voom_scale_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "scale",
    correccion = T
  )
print("######### SEXTA ########")

res1_voom_quantile_none_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "quantile",
    correccion = F
  )

res1_voom_quantile_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "quantile",
    correccion = T
  )
print("######### SEPTIMA ########")

res1_voom_quantile_none_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "cyclicloess",
    correccion = F
  )
print("######### OCTAVA ########")

res1_voom_none_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "cyclicloess",
    correccion = T
  )
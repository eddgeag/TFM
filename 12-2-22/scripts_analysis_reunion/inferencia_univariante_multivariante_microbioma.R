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
           correccion,
           microbiota.abs) {
    mc.instances <- numMCInstances(microbiota)
    mc.all <- getMonteCarloInstances(microbiota)
    res <- vector("list", length = mc.instances)
    res <- vector("list", length = mc.instances)
    
    vector.obesity <-
      vector("numeric", length = mc.instances)
    vector.grupo <-
      vector("numeric", length = mc.instances)
    vector.interaccion <-
      vector("numeric", length = mc.instances)
    
    
    for (mc.i in 1:mc.instances) {
      print(mc.i)
      t.input <- t(sapply(mc.all, function(y) {
        y[, mc.i]
      }))
      
      if (vom) {
        t.input <- t(t.input - min(t.input))
        t.input <- t(voom(t.input, normalize.method = "none")$E)
        
        
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
    
    
    pvaloresaux <-
      t(apply(Reduce("cbind", pvalores), 1, as.numeric))
    colnames(pvaloresaux) <-
      rep(colnames(pvalores[[1]]), mc.instances)
    diferenciasaux <-
      t(apply(Reduce("cbind", diferencias), 1, as.numeric))
    
    colnames(diferenciasaux) <-
      rep(colnames(diferencias[[1]]), mc.instances)
    
    matriz_pval <-
      matrix(NA,
             nrow = nrow(microbiota@reads),
             ncol = ncol(pvalores[[1]]))
    matriz_diferencias <-
      matrix(NA,
             nrow = nrow(microbiota@reads),
             ncol = ncol(diferencias[[1]]))
    
    nombres <- colnames(pvalores[[1]])
    nombresdif <- colnames(diferencias[[1]])
    
    for (col in 1:ncol(pvalores[[1]])) {
      nombre_i <- nombres[col]
      matriz_pval[, col] <-
        rowMeans(pvaloresaux[, grep(nombre_i, colnames(pvaloresaux))], na.rm = T)
      
    }
    rownames(matriz_pval) <- rownames(microbiota@reads)
    colnames(matriz_pval) <- colnames(pvalores[[1]])
    
    for (col in 1:ncol(diferencias[[1]])) {
      nombresdif_i <- nombresdif[col]
      matriz_diferencias[, col] <-
        rowMeans(diferenciasaux[, grep(nombresdif_i, x = colnames(diferenciasaux))], na.rm = T)
      
    }
    
    rownames(matriz_diferencias) <- rownames(microbiota@reads)
    colnames(matriz_diferencias) <- colnames(diferencias[[1]])
    adonis.obesidad <- mean(vector.obesity)
    adonis.grupo <- mean(vector.grupo)
    adonis.interaccion <- mean(vector.interaccion)
    
    matriz_prueba <-
      res[[1]][, grep("prueba", ignore.case = T, colnames(res[[1]]))]
    matriz_prueba <-
      matrix("MC_aldex",
             nrow = nrow(matriz_prueba),
             ncol = ncol(matriz_prueba))
    colnames(matriz_prueba) <-
      colnames(res[[1]])[grep("prueba", ignore.case = T, colnames(res[[1]]))]
    modelo <-
      bind_cols(matriz_prueba, matriz_diferencias, matriz_pval)
    resultado <- funcion_interpretacion(modelo, microbiota.abs)
    retorno <-
      list(
        resultado = resultado,
        adonis.interaccion = adonis.interaccion,
        adonis.grupo = adonis.grupo,
        adonis.obesidad = adonis.obesidad
      )
    return(retorno)
  }


microbiota <- datos$datos.filo
microbiota.abs <- datos$datos.filo.abs
res1_no_voom_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = F,
    method = "none",
    correccion = F,
    microbiota.abs = microbiota.abs
    
  )

print("######### PRIMERA ########")
res1_no_voom_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = F,
    method = "none",
    correccion = T,
    microbiota.abs = microbiota.abs
  )
print("######### SEGUNDA ########")

res1_voom_none_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "none",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### TERCERA ########")

res1_voom_none_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "none",
    correccion = T,
    microbiota.abs = microbiota.abs
  )

print("######### CUARTA ########")

res1_voom_scale_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "scale",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### QUINTA ########")

res1_voom_scale_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "scale",
    correccion = T,
    microbiota.abs = microbiota.abs
  )
print("######### SEXTA ########")

res1_voom_quantile_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "quantile",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### SEPTIMA ########")

res1_voom_quantile_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "quantile",
    correccion = T,
    microbiota.abs = microbiota.abs
  )
print("######### OCTAVA ########")

res1_voom_cyclicloess_none <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "cyclicloess",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### NOVENA ########")

res1_voom_cyclicloess_BH <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "cyclicloess",
    correccion = T,
    microbiota.abs = microbiota.abs
  )

microbiota <- datos$datos.genero
microbiota.abs <- datos$datos.genero.abs
res1_no_voom_none_genero <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = F,
    method = "none",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### PRIMERA ########")

res1_no_voom_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = F,
    method = "none",
    correccion = T,
    microbiota.abs = microbiota.abs
  )
print("######### SEGUNDA ########")

res1_voom_none_none_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "none",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### TERCERA ########")

res1_voom_none_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "none",
    correccion = T,
    microbiota.abs = microbiota.abs
  )

print("######### CUARTA ########")

res1_voom_scale_none_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "scale",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### QUINTA ########")

res1_voom_scale_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "scale",
    correccion = T,
    microbiota.abs = microbiota.abs
  )
print("######### SEXTA ########")

res1_voom_quantile_none_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "quantile",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### SEPTIMA ########")

res1_voom_quantile_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "quantile",
    correccion = T,
    microbiota.abs = microbiota.abs
  )
print("######### OCTAVA ########")

res1_voom_cyclicloess_none_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "cyclicloess",
    correccion = F,
    microbiota.abs = microbiota.abs
  )
print("######### NOVENA ########")

res1_voom_cyclicloess_BH_genero  <-
  analisis_microbiota(
    microbiota = microbiota,
    grupo = grupo,
    obesidad = obesidad,
    vom = T,
    method = "cyclicloess",
    correccion = T,
    microbiota.abs = microbiota.abs
  )



###### resumen y escritura de resultados #####
directorio_resultados <-
  file.path(ruta_directorio, "inferencia_univariante_metagenoma")
directorio_resultados_genero <-
  file.path(directorio_resultados, "genero")
directorio_resultados_BH_genero <-
  file.path(directorio_resultados_genero, "directorio_BH")
directorio_resultados_uncorrected_genero <-
  file.path(directorio_resultados_genero, "directorio_uncorrected")
directorio_resultados_filo <-
  file.path(directorio_resultados, "filo")
directorio_resultados_BH_filo <-
  file.path(directorio_resultados_filo, "directorio_BH")
directorio_resultados_uncorrected_filo <-
  file.path(directorio_resultados_filo, "directorio_uncorrected")

if (!dir.exists(directorio_resultados)) {
  dir.create(directorio_resultados)
}
if (!dir.exists(directorio_resultados_genero)) {
  dir.create(directorio_resultados_genero)
}
if (!dir.exists(directorio_resultados_BH_genero)) {
  dir.create(directorio_resultados_BH_genero)
}
if (!dir.exists(directorio_resultados_uncorrected_genero)) {
  dir.create(directorio_resultados_uncorrected_genero)
}
if (!dir.exists(directorio_resultados_filo)) {
  dir.create(directorio_resultados_filo)
}
if (!dir.exists(directorio_resultados_BH_filo)) {
  dir.create(directorio_resultados_BH_filo)
}
if (!dir.exists(directorio_resultados_uncorrected_filo)) {
  dir.create(directorio_resultados_uncorrected_filo)
}


print("######### analisis multivariante comparacion ##########")
print("FILO")

no_voom <- c(
  obesidad = res1_no_voom_BH$adonis.obesidad,
  grupo = res1_no_voom_BH$adonis.grupo,
  interaccion = res1_no_voom_BH$adonis.interaccion
)

voom_none <- c(
  obesidad = res1_voom_none_BH$adonis.obesidad,
  grupo = res1_voom_none_BH$adonis.grupo,
  interaccion = res1_voom_none_BH$adonis.interaccion
)

voom_scale <- c(
  obesidad = res1_voom_scale_BH$adonis.obesidad,
  grupo = res1_voom_scale_BH$adonis.grupo,
  interaccion = res1_voom_scale_BH$adonis.interaccion
)

voom_quantile <- c(
  obesidad = res1_voom_quantile_BH$adonis.obesidad,
  grupo = res1_voom_quantile_BH$adonis.grupo,
  interaccion = res1_voom_quantile_BH$adonis.interaccion
)

voom_cyclicloess <-
  c(
    obesidad = res1_voom_cyclicloess_BH$adonis.obesidad,
    grupo = res1_voom_cyclicloess_BH$adonis.grupo,
    interaccion = res1_voom_cyclicloess_BH$adonis.interaccion
  )


resultado_multivariante_filo <- rbind(no_voom,
                                      voom_none,
                                      voom_scale,
                                      voom_quantile,
                                      voom_cyclicloess)

write.csv(
  resultado_multivariante_filo,
  file = file.path(directorio_resultados_filo, "multivariante_filo.csv")
)

print("GENERO")

no_voom_genero <- c(
  obesidad = res1_no_voom_BH_genero$adonis.obesidad,
  grupo = res1_no_voom_BH_genero$adonis.grupo,
  interaccion = res1_no_voom_BH_genero$adonis.interaccion
)

voom_none_genero <-
  c(
    obesidad = res1_voom_none_BH_genero$adonis.obesidad,
    grupo = res1_voom_none_BH_genero$adonis.grupo,
    interaccion = res1_voom_none_BH_genero$adonis.interaccion
  )

voom_scale_genero <-
  c(
    obesidad = res1_voom_scale_BH_genero$adonis.obesidad,
    grupo = res1_voom_scale_BH_genero$adonis.grupo,
    interaccion = res1_voom_scale_BH_genero$adonis.interaccion
  )

voom_quantile_genero <-
  c(
    obesidad = res1_voom_quantile_BH_genero$adonis.obesidad,
    grupo = res1_voom_quantile_BH_genero$adonis.grupo,
    interaccion = res1_voom_quantile_BH_genero$adonis.interaccion
  )

voom_cyclicloess_genero <-
  c(
    obesidad = res1_voom_cyclicloess_BH_genero$adonis.obesidad,
    grupo = res1_voom_cyclicloess_BH_genero$adonis.grupo,
    interaccion = res1_voom_cyclicloess_BH_genero$adonis.interaccion
  )

resultado_multivariante_genero <- rbind(
  no_voom_genero,
  voom_none_genero,
  voom_scale_genero,
  voom_quantile_genero,
  voom_cyclicloess_genero
)
write.csv(
  resultado_multivariante_genero,
  file = file.path(directorio_resultados_genero, "multivariante_genero.csv")
)

print("######### ESCRITURA RESULTADOS FILO ############")
## no voom
write.csv(
  res1_no_voom_none$resultado,
  file.path(
    directorio_resultados_uncorrected_filo,
    "filo_no_voom_none.csv"
  )
)

write.csv(
  res1_no_voom_BH$resultado,
  file.path(directorio_resultados_BH_filo,
            "filo_no_voom_none.csv")
)
## voom no method
write.csv(
  res1_voom_none_none$resultado,
  file.path(
    directorio_resultados_uncorrected_filo,
    "filo_voom_none.csv"
  )
)

write.csv(
  res1_voom_none_BH$resultado,
  file.path(directorio_resultados_BH_filo,
            "filo_voom_none.csv")
)
## voom scaled
write.csv(
  res1_voom_scale_none$resultado,
  file.path(
    directorio_resultados_uncorrected_filo,
    "filo_voom_scaled.csv"
  )
)

write.csv(
  res1_voom_scale_BH$resultado,
  file.path(directorio_resultados_BH_filo,
            "filo_voom_scaled.csv")
)
## voom quantile
write.csv(
  res1_voom_quantile_none$resultado,
  file.path(
    directorio_resultados_uncorrected_filo,
    "filo_voom_quantile.csv"
  )
)
write.csv(
  res1_voom_quantile_BH$resultado,
  file.path(directorio_resultados_BH_filo,
            "filo_voom_quantile.csv")
)

## voom cycloless
write.csv(
  res1_voom_cyclicloess_none$resultado,
  file.path(
    directorio_resultados_uncorrected_filo,
    "filo_voom_cyclicloess.csv"
  )
)
write.csv(
  res1_voom_cyclicloess_BH$resultado,
  file.path(directorio_resultados_BH_filo,
            "filo_voom_cyclicloess.csv")
)


print("######### ESCRITURA RESULTADOS GENERO ############")
## no voom
write.csv(
  res1_no_voom_none_genero$resultado,
  file.path(
    directorio_resultados_uncorrected_genero,
    "filo_no_voom_none.csv"
  )
)

write.csv(
  res1_no_voom_BH_genero$resultado,
  file.path(directorio_resultados_BH_genero,
            "filo_no_voom_none.csv")
)
## voom no method
write.csv(
  res1_voom_none_none_genero$resultado,
  file.path(
    directorio_resultados_uncorrected_genero,
    "filo_voom_none.csv"
  )
)

write.csv(
  res1_voom_none_BH_genero$resultado,
  file.path(directorio_resultados_BH_genero,
            "filo_voom_none.csv")
)
## voom scaled
write.csv(
  res1_voom_scale_none_genero$resultado,
  file.path(
    directorio_resultados_uncorrected_genero,
    "filo_voom_scaled.csv"
  )
)

write.csv(
  res1_voom_scale_BH_genero$resultado,
  file.path(directorio_resultados_BH_genero,
            "filo_voom_scaled.csv")
)
## voom quantile
write.csv(
  res1_voom_quantile_none_genero$resultado,
  file.path(
    directorio_resultados_uncorrected_genero,
    "filo_voom_quantile.csv"
  )
)
write.csv(
  res1_voom_quantile_BH_genero$resultado,
  file.path(directorio_resultados_BH_genero,
            "filo_voom_quantile.csv")
)

## voom cycloless
write.csv(
  res1_voom_cyclicloess_none_genero$resultado,
  file.path(
    directorio_resultados_uncorrected_genero,
    "filo_voom_cyclicloess.csv"
  )
)
write.csv(
  res1_voom_cyclicloess_BH_genero$resultado,
  file.path(
    directorio_resultados_BH_genero,
    "filo_voom_cyclicloess.csv"
  )
)

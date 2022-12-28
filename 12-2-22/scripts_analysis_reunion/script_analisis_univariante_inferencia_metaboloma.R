
# EN ESTE SCRIPT ANALIZAREMOS EL METABOLOMA

# LLAMAMOS AL SCRIPT DE LA INFERENCIA UNIVARIANTE

add_classificacion_metaboloma <- function(X){
  
  aromatic <- colnames(metaboloma)[c(2, 3, 4, 1, 6, 35, 36, 34)]
  other <- colnames(metaboloma)[c(21, 9, 16, 26, 22, 18, 13, 29, 30, 27)]
  aa <- colnames(metaboloma)[c(33, 15, 20, 19, 23, 25)]
  carbo <-
    colnames(metaboloma)[!colnames(metaboloma) %in% c(aromatic, other, aa)]
  
  X$grupo_meta <-
    ifelse(colnames(metaboloma) %in% aromatic, "AROMATICO", NA)
  X$grupo_meta <-
    ifelse(colnames(metaboloma) %in% other,
           "OTHER",
           X$grupo_meta)
  X$grupo_meta <-
    ifelse(colnames(metaboloma) %in% aa,
           "DERIVADOS AA",
           X$grupo_meta)
  X$grupo_meta <-
    ifelse(
      colnames(metaboloma) %in% carbo,
      "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
      X$grupo_meta
    )
  
  return(X)
}

construct_pre_Results <- function(X,correct,grupo,obesidad){
  
  resultado <- inferencia_univariante(X,grupo = grupo,obesidad=obesidad)
  p.values <- resultado$p.values
  diferencias <- resultado$diferencias
  tipos_de_pruebas <- resultado$tipo_prueba
  if(!correct){
    
  }else{
    
    p.values[,1:4] <- apply(p.values[,1:5],2,function(x) p.adjust(x,"BH"))
    p.values[,5:ncol(p.values)] <- apply(p.values[,5:ncol(p.values)],1,function(x) p.adjust(x,"BH"))
    
  }
  modelo <- bind_cols(tipos_de_pruebas,diferencias,p.values)
  
  return(modelo)
  
}

funcion_wrapper <- function(X,grupo_,obesidad_){
  
  modelo <- construct_pre_Results(X,correct = F,grupo = grupo_,obesidad = obesidad_)
  modelo.BH <- construct_pre_Results(X,correct = T,grupo = grupo_,obesidad = obesidad_)
  
  res1 <-add_classificacion_metaboloma(funcion_interpretacion(glm_modelo= modelo,X))
  res1 <- res1[order(res1$grupo_meta),]
  res2 <-add_classificacion_metaboloma(funcion_interpretacion(glm_modelo= modelo.BH,X))
  res2 <- res2[order(res2$grupo_meta),]
  
  return(list(correctedBH=res2,normal=res1))
}
# LEEMOS LOS DATOS

ruta_directorio <- "../../datos/resultados_analisis_datos_procesados"
ruta_archivo_datos <- "datos_integromics_procesados.rds"
ruta_actual_directorio <- "./12-2-22/scripts_analysis_reunion"
funciones_inferencia <- "infererencia_univariante_funcion.R"
datos <- readRDS(file.path(ruta_directorio,ruta_archivo_datos))
ruta_directorio_resultados_metaboloma <- "inferencia_univariante_metaboloma"
directorio_resultados <- file.path(ruta_directorio,
                                   ruta_directorio_resultados_metaboloma)
if(!dir.exists(directorio_resultados)){
  dir.create(directorio_resultados)
}
directorio_sin_corregir <- file.path(directorio_resultados,"not_corrected")
directorio_corregido <- file.path(directorio_resultados,"BH_corrected")

if(!dir.exists(directorio_sin_corregir)){
  dir.create(directorio_sin_corregir)
}
if(!dir.exists(directorio_corregido)){
  dir.create(directorio_corregido)
}

#  necesitamos general data y variables_in_bacteria
general_data <- datos$datos.clinicos.totales
variables_in_bacteria <- datos$variables_in_bacteria
# necesitamos las variables dependientes
grupo <- datos$grupo.total
obesidad <- datos$obesidad.total

source(file.path(ruta_actual_directorio,funciones_inferencia))


### metaboloma 53 no voom
metaboloma <- datos$datos.metaboloma.totales
grupo <- datos$grupo.total
obesidad <- datos$obesidad.total
resultado_53_no_voom <- funcion_wrapper(metaboloma,grupo,obesidad)
write.csv(apply(resultado_53_no_voom$normal,2,as.character),file=file.path(directorio_sin_corregir,"53_no_voom_none.csv"),quote = F,row.names = F)
write.csv(apply(resultado_53_no_voom$correctedBH,2,as.character),file=file.path(directorio_corregido,"53_no_voom_BH.csv"),quote = F,row.names = F)

### metaboloma 53 voom

metaboloma <- datos$datos.metaboloma.totales.voom
grupo <- datos$grupo.total
obesidad <- datos$obesidad.total
resultado_53_voom <- funcion_wrapper(metaboloma,grupo,obesidad)
write.csv(apply(resultado_53_voom$normal,2,as.character),file=file.path(directorio_sin_corregir,"53_voom_none.csv"),quote = F,row.names = F)
write.csv(apply(resultado_53_voom$correctedBH,2,as.character),file=file.path(directorio_corregido,"53_voom_BH.csv"),quote = F,row.names = F)

### metaboloma 46 common no voom

metaboloma <- datos$datos.metaboloma.common
grupo <- datos$grupo.common
obesidad <- datos$obesidad.common
resultado_46_no_voom <- funcion_wrapper(metaboloma,grupo,obesidad)
write.csv(apply(resultado_46_no_voom$normal,2,as.character),file=file.path(directorio_sin_corregir,"46_no_voom_none.csv"),quote = F,row.names = F)
write.csv(apply(resultado_46_no_voom$correctedBH,2,as.character),file=file.path(directorio_corregido,"46_no_voom_BH.csv"),quote = F,row.names = F)

### metaboloma 46 common voom

metaboloma <- datos$datos.metaboloma.common.voom
grupo <- datos$grupo.common
obesidad <- datos$obesidad.common
resultado_46_voom <- funcion_wrapper(metaboloma,grupo,obesidad)
write.csv(apply(resultado_46_voom$normal,2,as.character),file=file.path(directorio_sin_corregir,"46_voom_none.csv"),quote = F,row.names = F)
write.csv(apply(resultado_46_voom$correctedBH,2,as.character),file=file.path(directorio_corregido,"46_voom_BH.csv"),quote = F,row.names = F)






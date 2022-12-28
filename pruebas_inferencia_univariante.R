




fn.paired <- function(x, g, termino1, termino2) {
  x1 <- x[g == termino1]
  x2 <- x[g == termino2]
  mod <- lm(x ~ g)
  prueba.normalidad <- shapiro.test(residuals(mod))$p.value
  prueba.homocedasticidad <- leveneTest(mod)$`Pr(>F)`[1]
  if (length(x1) == 0 || length(x2) == 0) {
    if (length(x1) == 0) {
      p.prueba <- 0
      p.location <- mean(x1)
      tipo.prueba <- paste0("Ausencia de ", termino1)
      
    } else if (length(x2) == 0) {
      p.prueba <- 0
      p.location <- mean(x2)
      tipo.prueba < -paste0("Ausencia de ", termino2)
      
    }
  }
  else if (prueba.normalidad < 0.05) {
    p.prueba <- wilcox.test(x ~ g)$p.value
    p.location <- diff(c(median(x1), median(x2)))
    tipo.prueba <- "noparam"
    
  } else if (prueba.normalidad > 0.05 &&
             prueba.homocedasticidad > 0.05) {
    p.prueba <- t.test(x ~ g, var.equal = T)$p.value
    p.location <- diff(c(mean(x1), mean(x2)))
    p.prueba <- "param"
    
  } else if (prueba.normalidad > 0.05 &&
             prueba.homocedasticidad < 0.05) {
    p.prueba <- t.test(x ~ g, var.equal = F)$p.value
    p.location <- diff(c(mean(x1), mean(x2)))
    tipo.prueba <- "param"
  }
  
  return(p.prueba = p.prueba,
         p.location = p.location,
         tipo.prueba = tipo.prueba)
}


X <- metaboloma

p <- ncol(X)
n <- nrow(X)

if (n == 46) {
  grupo <- variables_in_bacteria$GROUP
  obesidad <- variables_in_bacteria$OBESE
} else{
  grupo <- general_data$GROUP
  obesidad <- general_data$OBESE
}
### VARIABLES
obesidad.mod <- vector("numeric", length = n)
s.obesidad.mod1 <- vector("numeric", length = n)
s.grupo.mod1 <- vector("numeric", length = n)
s.interaccion.mod1 <- vector("numeric", length = n)
s.prueba.mod1 <-  vector("character", length = n)

s.obesidad.mod2 <- vector("numeric", length = n)
s.grupo.mod2 <- vector("numeric", length = n)
s.interaccion.mod2 <- vector("numeric", length = n)
s.prueba.mod2 <-  vector("character", length = n)

obesidad.pvalue <-  vector("numeric", length = n)
diferencia_loc_obesidad <- vector("numeric", length = n)
obesidad.prueba <- vector("character", length = n)

obesidad.pvalue2 <-  vector("numeric", length = n)
diferencia_loc_obesidad2 <- vector("numeric", length = n)
obesidad.prueba2 <- vector("character", length = n)

n.sexo <- length(grupo[grupo != "PCOS"])
sexo.pvalue <-  vector("numeric", length = n.sexo)
diferencia_loc_sexo <- vector("numeric", length = n.sexo)
sexo.prueba <- vector("character", length = n.sexo)

n.malespcos <- length(grupo[grupo != "Female"])
malespcos.pvalue <-  vector("numeric", length = n.malespcos)
diferencia_loc_malespcos <- vector("numeric", length = n.malespcos)
malespcos.prueba <- vector("character", length = n.malespcos)

n.femalespcos <- length(grupo[grupo != "Male"])
femalespcos.pvalue <-  vector("numeric", length = n.femalespcos)
diferencia_loc_femalespcos <-
  vector("numeric", length = n.femalespcos)
femalespcos.prueba <- vector("character", length = n.femalespcos)

sexo.pvalue2 <-  vector("numeric", length = n.sexo)
diferencia_loc_sexo2 <- vector("numeric", length = n.sexo)
sexo.prueba2 <- vector("character", length = n.sexo)


malespcos.pvalue2 <-  vector("numeric", length = n.malespcos)
diferencia_loc_malespcos2 <- vector("numeric", length = n.malespcos)
malespcos.prueba2 <- vector("character", length = n.malespcos)


femalespcos.pvalue2 <-  vector("numeric", length = n.femalespcos)
diferencia_loc_femalespcos2 <-
  vector("numeric", length = n.femalespcos)
femalespcos.prueba2 <- vector("character", length = n.femalespcos)


n.sexo.no.obese <-
  length(grupo[grupo != "PCOS" & obesidad == "No Obese"])
sexo.pvalue.no.obese <-  vector("numeric", length = n.sexo.no.obese)
diferencia_loc_sexo.no.obese <-
  vector("numeric", length = n.sexo.no.obese)
sexo.prueba.no.obese <-
  vector("character", length = n.sexo.no.obese)

n.malespcos.no.obese <-
  length(grupo[grupo != "Female" & obesidad == "No Obese"])
malespcos.pvalue.no.obese <-
  vector("numeric", length = n.malespcos.no.obese)
diferencia_loc_malespcos.no.obese <-
  vector("numeric", length = n.malespcos.no.obese)
malespcos.prueba.no.obese <-
  vector("character", length = n.malespcos.no.obese)

n.femalespcos.no.obese <-
  length(grupo[grupo != "Male" & obesidad == "No Obese"])
femalespcos.pvalue.no.obese <-
  vector("numeric", length = n.femalespcos.no.obese)
diferencia_loc_femalespcos.no.obese <-
  vector("numeric", length = n.femalespcos.no.obese)
femalespcos.prueba.no.obese <-
  vector("character", length = n.femalespcos.no.obese)


n.sexo.obese <- length(grupo[grupo != "PCOS" & obesidad == "Obese"])
sexo.pvalue.obese <-  vector("numeric", length = n.sexo.obese)
diferencia_loc_sexo.obese <-
  vector("numeric", length = n.sexo.obese)
sexo.prueba.obese <- vector("character", length = n.sexo.obese)

n.malespcos.obese <-
  length(grupo[grupo != "Female" & obesidad == "Obese"])
malespcos.pvalue.obese <-
  vector("numeric", length = n.malespcos.obese)
diferencia_loc_malespcos.obese <-
  vector("numeric", length = n.malespcos.obese)
malespcos.prueba.obese <-
  vector("character", length = n.malespcos.obese)

n.femalespcos.obese <-
  length(grupo[grupo != "Male" & obesidad == "Obese"])
femalespcos.pvalue.obese <-
  vector("numeric", length = n.femalespcos.obese)
diferencia_loc_femalespcos.obese <-
  vector("numeric", length = n.femalespcos.obese)
femalespcos.prueba.obese <-
  vector("character", length = n.femalespcos.obese)



n.sexo.no.obese <-
  length(grupo[grupo != "PCOS" & obesidad == "No Obese"])
sexo.pvalue.no.obese2 <-
  vector("numeric", length = n.sexo.no.obese)
diferencia_loc_sexo.no.obese2 <-
  vector("numeric", length = n.sexo.no.obese)
sexo.prueba.no.obese2 <-
  vector("character", length = n.sexo.no.obese)

n.malespcos.no.obese <-
  length(grupo[grupo != "Female" & obesidad == "No Obese"])
malespcos.pvalue.no.obese2 <-
  vector("numeric", length = n.malespcos.no.obese)
diferencia_loc_malespcos.no.obese2 <-
  vector("numeric", length = n.malespcos.no.obese)
malespcos.prueba.no.obese2 <-
  vector("character", length = n.malespcos.no.obese)

n.femalespcos.no.obese <-
  length(grupo[grupo != "Male" & obesidad == "No Obese"])
femalespcos.pvalue.no.obese2 <-
  vector("numeric", length = n.femalespcos.no.obese)
diferencia_loc_femalespcos.no.obese2 <-
  vector("numeric", length = n.femalespcos.no.obese)
femalespcos.prueba.no.obese2 <-
  vector("character", length = n.femalespcos.no.obese)


n.sexo.obese <- length(grupo[grupo != "PCOS" & obesidad == "Obese"])
sexo.pvalue.obese2 <-  vector("numeric", length = n.sexo.obese)
diferencia_loc_sexo.obese2 <-
  vector("numeric", length = n.sexo.obese)
sexo.prueba.obese2 <- vector("character", length = n.sexo.obese)

n.malespcos.obese <-
  length(grupo[grupo != "Female" & obesidad == "Obese"])
malespcos.pvalue.obese2 <-
  vector("numeric", length = n.malespcos.obese)
diferencia_loc_malespcos.obese2 <-
  vector("numeric", length = n.malespcos.obese)
malespcos.prueba.obese2 <-
  vector("character", length = n.malespcos.obese)

n.femalespcos.obese <-
  length(grupo[grupo != "Male" & obesidad == "Obese"])
femalespcos.pvalue.obese2 <-
  vector("numeric", length = n.femalespcos.obese)
diferencia_loc_femalespcos.obese2 <-
  vector("numeric", length = n.femalespcos.obese)
femalespcos.prueba.obese2 <-
  vector("character", length = n.femalespcos.obese)



n.sexo.obesidad <- length(obesidad[grupo != "PCOS"])
sexo.p.value.obesidad <-
  vector("numeric", length = n.sexo.obesidad)
sexo.diferencia.obesidad <-
  vector("numeric", length = n.sexo.obesidad)
sexo.prueba.obesidad <-
  vector("character", length = n.sexo.obesidad)

n.fems.obesidad <- length(obesidad[grupo != "Male"])
fems.p.value.obesidad <- vector("numeric", length = n.fems.obesidad)
fems.diferencia.obesidad <-
  vector("numeric", length = n.fems.obesidad)
fems.prueba.obesidad <- vector("numeric", length = n.fems.obesidad)

n.pcosmales.obesidad <- length(obesidad[grupo != "Female"])
pcosmales.p.value.obesidad <-
  vector("numeric", length = n.pcosmales.obesidad)
pcosmales.diferencia.obesidad <-
  vector("numeric", length = n.pcosmales.obesidad)
pcosmales.prueba.obesidad <-
  vector("numeric", length = n.pcosmales.obesidad)


n.sexo.obesidad <- length(obesidad[grupo != "PCOS"])
sexo.p.value.obesidad2 <-
  vector("numeric", length = n.sexo.obesidad)
sexo.diferencia.obesidad2 <-
  vector("numeric", length = n.sexo.obesidad)
sexo.prueba.obesidad2 <-
  vector("character", length = n.sexo.obesidad)

n.fems.obesidad <- length(obesidad[grupo != "Male"])
fems.p.value.obesidad2 <-
  vector("numeric", length = n.fems.obesidad)
fems.diferencia.obesidad2 <-
  vector("numeric", length = n.fems.obesidad)
fems.prueba.obesidad2 <- vector("numeric", length = n.fems.obesidad)

n.pcosmales.obesidad <- length(obesidad[grupo != "Female"])
pcosmales.p.value.obesidad2 <-
  vector("numeric", length = n.pcosmales.obesidad)
pcosmales.diferencia.obesidad2 <-
  vector("numeric", length = n.pcosmales.obesidad)
pcosmales.prueba.obesidad2 <-
  vector("numeric", length = n.pcosmales.obesidad)


for (j in 1:ncol(p)) {
  var_i <- X[, j]
  modelo1 <-
    lm(var_i ~ grupo * obesidad,
       contrasts = list(grupo = contr.sum, obesidad = contr.sum))
  modelo2 <-
    lm(var_i ~ obesidad * grupo,
       contrasts = list(grupo = contr.sum, obesidad = contr.sum))
  
  ### vemos las suposiciones de los modelos
  
  prueba.normalidad <- shapiro.test(residuals(modelo1))$p.value
  prueba.homocedasticidad <- leveneTest(modelo1)$`Pr(>F)`[1]
  prueba.normalidad2 <- shapiro.test(residuals(modelo2))$p.value
  prueba.homocedasticidad2 <- leveneTest(modelo2)$`Pr(>F)`[1]
  ### si es cumple los criterios de normalidad
  
  if (prueba.normalidad < 0.05) {
    s <-
      summary(aovp(
        var_i ~ grupo * obesidad,
        contrasts = list(grupo = contr.sum, obesidad = contr.sum)
      ))
    s.obesidad.mod1[j] <- s[[1]]$`Pr(Prob)`[2]
    s.grupo.mod1[j] <- s[[1]]$`Pr(Prob)`[1]
    s.interaccion.mod1[j] <- s[[1]]$`Pr(Prob)`[3]
    s.prueba.mod1[j] <- "noparam"
    
    
  } else{
    s <- Anova(modelo1, type = 3)
    s.obesidad.mod1[j] <- s[[1]]$`Pr(Prob)`[2]
    s.grupo.mod1[j] <- s[[1]]$`Pr(Prob)`[1]
    s.interaccion.mod1[j] <- s[[1]]$`Pr(Prob)`[3]
    s.prueba.mod1[j] <- "param"
    
  }
  if (prueba.normalidad2 < 0.05) {
    s <-
      summary(aovp(
        var_i ~ obesidad * grupo,
        contrasts = list(grupo = contr.sum, obesidad = contr.sum)
      ))
    s.obesidad.mod2[j] <- s[[1]]$`Pr(Prob)`[2]
    s.grupo.mod2[j] <- s[[1]]$`Pr(Prob)`[1]
    s.interaccion.mod2[j] <- s[[1]]$`Pr(Prob)`[3]
    s.prueba.mod2[j] <- "noparam"
    
    
  } else{
    s <- Anova(modelo2, type = 3)
    s.obesidad.mod2[j] <- s[[1]]$`Pr(Prob)`[2]
    s.grupo.mod2[j] <- s[[1]]$`Pr(Prob)`[1]
    s.interaccion.mod2[j] <- s[[1]]$`Pr(Prob)`[3]
    s.prueba.mod2[j] <- "param"
    
  }
  
  ### veamos la obesidad del primer modelo
  
  if (s.obesidad.mod1[j] < 0.05) {
    vector_resultados <-
      fn.paired(
        x = var_i,
        g = obesidad,
        termino1 = "No Obese",
        termino2 = "Obese"
      )
    obesidad.pvalue[j] = vector_resultados[1]
    diferencia_loc_obesidad[j] = vector_resultados[2]
    obesidad.prueba[j] = vector_resultados[3]
  }
  
  ### veamos la obesidad del segundo modelo
  
  if (s.obesidad.mod2[j] < 0.05) {
    vector_resultados <-
      fn.paired(
        x = var_i,
        g = obesidad,
        termino1 = "No Obese",
        termino2 = "Obese"
      )
    obesidad.pvalue2[j] = vector_resultados[1]
    diferencia_loc_obesidad2[j] = vector_resultados[2]
    obesidad.prueba2[j] = vector_resultados[3]
    
    
  }
  ### grupo modelo 1
  if (s.grupo.mod1[j] < 0.05) {
    ### sexo
    var_i_sexo <- var_i[grupo != "PCOS"]
    grupo_sexo <- grupo[grupo != "PCOS"]
    vector_resultados <-
      fn.paired(
        x = var_i_sexo,
        g = grupo_sexo,
        termino1 = "Female",
        termino2 = "Male"
      )
    sexo.pvalue[j] = vector_resultados[1]
    diferencia_loc_sexo[j] = vector_resultados[2]
    sexo.prueba[j] = vector_resultados[3]
    ### hombres vs mujeres PCOS
    
    var_i_malespcos <- var_i[grupo != "Female"]
    grupo_malespcos <- grupo[grupo != "Female"]
    vector_resultados <-
      fn.paired(
        x = var_i_malespcos,
        g = grupo_sexo,
        termino1 = "Male",
        termino2 = "PCOS"
      )
    malespcos.pvalue[j] = vector_resultados[1]
    diferencia_loc_malespcos[j] = vector_resultados[2]
    malespcos.prueba[j] = vector_resultados[3]
    
    ### mujeres vs mujeres PCOS
    var_i_femalespcos <- var_i[grupo != "Male"]
    grupo_femalespcos <- grupo[grupo != "Male"]
    vector_resultados <-
      fn.paired(
        x = var_i_femalespcos,
        g = grupo_sexo,
        termino1 = "Female",
        termino2 = "PCOS"
      )
    femalespcos.pvalue[j] = vector_resultados[1]
    diferencia_loc_femalespcos[j] = vector_resultados[2]
    femalespcos.prueba[j] = vector_resultados[3]
    
  }
  #### grupo modelo 2
  if (s.grupo.mod2[j] < 0.05) {
    ### sexo
    var_i_sexo <- var_i[grupo != "PCOS"]
    grupo_sexo <- grupo[grupo != "PCOS"]
    vector_resultados <-
      fn.paired(
        x = var_i_sexo,
        g = grupo_sexo,
        termino1 = "Female",
        termino2 = "Male"
      )
    sexo.pvalue2[j] = vector_resultados[1]
    diferencia_loc_sexo2[j] = vector_resultados[2]
    sexo.prueba2[j] = vector_resultados[3]
    ### hombres vs mujeres PCOS
    
    var_i_malespcos <- var_i[grupo != "Female"]
    grupo_malespcos <- grupo[grupo != "Female"]
    vector_resultados <-
      fn.paired(
        x = var_i_malespcos,
        g = grupo_sexo,
        termino1 = "Male",
        termino2 = "PCOS"
      )
    malespcos.pvalue2[j] = vector_resultados[1]
    diferencia_loc_malespcos2[j] = vector_resultados[2]
    malespcos.prueba2[j] = vector_resultados[3]
    
    ### mujeres vs mujeres PCOS
    var_i_femalespcos <- var_i[grupo != "Male"]
    grupo_femalespcos <- grupo[grupo != "Male"]
    vector_resultados <-
      fn.paired(
        x = var_i_femalespcos,
        g = grupo_sexo,
        termino1 = "Female",
        termino2 = "PCOS"
      )
    femalespcos.pvalue2[j] = vector_resultados[1]
    diferencia_loc_femalespcos2[j] = vector_resultados[2]
    femalespcos.prueba2[j] = vector_resultados[3]
    
  }
  if (s.interaccion.mod1[j] < 0.05) {
    ### obesidad|grupo
    ### Vemos el sexo en las personas no obesas
    
    No.Obese_var_i_sexo <-
      var_i[obesidad == "No Obese" & grupo != "PCOS"]
    No.Obese_grupo_sexo <-
      grupo[obesidad == "No Obese" & grupo != "PCOS"]
    vector_resultados <-
      fn.paired(
        x = No.Obese_var_i_sexo,
        g = No.Obese_grupo_sexo,
        termino1 = "Female",
        termino2 = "Male"
      )
    sexo.pvalue.no.obese[j] = vector_resultados[1]
    diferencia_loc_sexono.obese[j] = vector_resultados[2]
    sexo.prueba.no.obese[j] = vector_resultados[3]
    ###### Vemos hombres vs mujeres PCOS no obesas
    
    No.Obese_var_i_malespcos <-
      var_i[obesidad == "No Obese" & grupo != "Female"]
    No.Obese_grupo_malespcos <-
      grupo[obesidad == "No Obese" & grupo != "Female"]
    vector_resultados <-
      fn.paired(
        x = No.Obese_var_i_malespcos,
        g = No.Obese_grupo_malespcos,
        termino1 = "Male",
        termino2 = "PCOS"
      )
    malespcos.pvalue.no.obese[j] = vector_resultados[1]
    diferencia_loc_malespcos.no.obese[j] = vector_resultados[2]
    malespcos.prueba.no.obese[j] = vector_resultados[3]
    
    ### No obesidad en mujeres vs mujeres PCOS
    No.Obese_var_i_femalespcos <-
      var_i[obesidad == "No Obese" & grupo != "Male"]
    No.Obese_grupo_femalespcos <-
      grupo[obesidad == "No Obese" & grupo != "Male"]
    vector_resultados <-
      fn.paired(
        x = No.Obese_var_i_femalespcos,
        g = No.Obese_grupo_femalespcos,
        termino1 = "Female",
        termino2 = "PCOS"
      )
    femalespcos.pvalue.no.obese[j] = vector_resultados[1]
    diferencia_loc_femalespcos.no.obese[j] = vector_resultados[2]
    femalespcos.prueba.no.obese[j] = vector_resultados[3]
    
    
    ### obesidad|grupo
    ### Vemos el sexo en las personas obesas
    
    Obese_var_i_sexo <- var_i[obesidad == "Obese" & grupo != "PCOS"]
    Obese_grupo_sexo <-  grupo[obesidad == "Obese" & grupo != "PCOS"]
    vector_resultados <-
      fn.paired(
        x = Obese_var_i_sexo,
        g = Obese_grupo_sexo,
        termino1 = "Female",
        termino2 = "Male"
      )
    sexo.pvalue.obese[j] = vector_resultados[1]
    diferencia_loc_sexo.obese[j] = vector_resultados[2]
    sexo.prueba.obese[j] = vector_resultados[3]
    ###### Vemos hombres vs mujeres PCOS obesas
    
    Obese_var_i_malespcos <-
      var_i[obesidad == "Obese" & grupo != "Female"]
    Obese_grupo_malespcos <-
      grupo[obesidad == "Obese" & grupo != "Female"]
    vector_resultados <-
      fn.paired(
        x = Obese_var_i_malespcos,
        g = Obese_grupo_malespcos,
        termino1 = "Male",
        termino2 = "PCOS"
      )
    malespcos.pvalue.obese[j] = vector_resultados[1]
    diferencia_loc_malespcos.obese[j] = vector_resultados[2]
    malespcos.prueba.obese[j] = vector_resultados[3]
    
    ### obesidad en mujeres vs mujeres PCOS
    Obese_var_i_femalespcos <-
      var_i[obesidad == "Obese" & grupo != "Male"]
    Obese_grupo_femalespcos <-
      grupo[obesidad == "Obese" & grupo != "Male"]
    vector_resultados <-
      fn.paired(
        x = Obese_var_i_femalespcos,
        g = Obese_grupo_femalespcos,
        termino1 = "Female",
        termino2 = "PCOS"
      )
    femalespcos.pvalue.obese[j] = vector_resultados[1]
    diferencia_loc_femalespcos.obese[j] = vector_resultados[2]
    femalespcos.prueba.obese[j] = vector_resultados[3]
    
    
    ### AHORA VEREMOS EL GRUPO|OBESIDAD
    
    sexo_var_i <- var_i[grupo != "PCOS"]
    sexo.obesidad <- obesidad[grupo != "PCOS"]
    
    fems_var_i <- var_i[grupo != "Male"]
    fems.obesidad <- obesidad[grupo != "Male"]
    
    pcosmales_var_i <- var_i[grupo != "Female"]
    pcosmales.obesidad <- obesidad[grupo != "Female"]
    
    vector_resultados <-
      fn.paired(
        x = sexo_var_i,
        g = sexo.obesidad,
        termino1 = "No Obese",
        termino2 = "Obese"
      )
    sexo.p.value.obesidad[j] <- vector_resultados[1]
    sexo.diferencia.obesidad[j]  <- vector_resultados[2]
    sexo.prueba.obesidad[j]  <- vector_resultados[3]
    
    vector_resultados <-
      fn.paired(
        x = fems_var_i,
        g = fems.obesidad,
        termino1 = "No Obese",
        termino2 = "Obese"
      )
    fems.p.value.obesidad[j]  <- vector_resultados[1]
    fems.diferencia.obesidad[j]  <- vector_resultados[2]
    fems.prueba.obesidad[j]  <- vector_resultados[3]
    
    
    vector_resultados <-
      fn.paired(
        x = pcosmales_var_i,
        g = pcosmales.obesidad,
        termino1 = "No Obese",
        termino2 = "Obese"
      )
    pcosmales.p.value.obesidad[j]  <- vector_resultados[1]
    pcosmales.diferencia.obesidad[j]  <- vector_resultados[2]
    pcosmales.prueba.obesidad[j]  <- vector_resultados[3]
    
    
    
  }
  
  if (s.interaccion.mod2[j] < 0.05) {
    ### obesidad|grupo
    ### Vemos el sexo en las personas no obesas
    
    No.Obese_var_i_sexo <-
      var_i[obesidad == "No Obese" & grupo != "PCOS"]
    No.Obese_grupo_sexo <-
      grupo[obesidad == "No Obese" & grupo != "PCOS"]
    vector_resultados <-
      fn.paired(
        x = No.Obese_var_i_sexo,
        g = No.Obese_grupo_sexo,
        termino1 = "Female",
        termino2 = "Male"
      )
    sexo.pvalue.no.obese2[j] = vector_resultados[1]
    diferencia_loc_sexo.no.obese2[j] = vector_resultados[2]
    sexo.prueba.no.obese2[j] = vector_resultados[3]
    ###### Vemos hombres vs mujeres PCOS no obesas
    
    No.Obese_var_i_malespcos <-
      var_i[obesidad == "No Obese" & grupo != "Female"]
    No.Obese_grupo_malespcos <-
      grupo[obesidad == "No Obese" & grupo != "Female"]
    vector_resultados <-
      fn.paired(
        x = No.Obese_var_i_malespcos,
        g = No.Obese_grupo_malespcos,
        termino1 = "Male",
        termino2 = "PCOS"
      )
    malespcos.pvalue.no.obese2[j] = vector_resultados[1]
    diferencia_loc_malespcos.no.obese2[j] = vector_resultados[2]
    malespcos.prueba.no.obese2[j] = vector_resultados[3]
    
    ### No obesidad en mujeres vs mujeres PCOS
    No.Obese_var_i_femalespcos <-
      var_i[obesidad == "No Obese" & grupo != "Male"]
    No.Obese_grupo_femalespcos <-
      grupo[obesidad == "No Obese" & grupo != "Male"]
    vector_resultados <-
      fn.paired(
        x = No.Obese_var_i_femalespcos,
        g = No.Obese_grupo_femalespcos,
        termino1 = "Female",
        termino2 = "PCOS"
      )
    femalespcos.pvalue.no.obese2[j] = vector_resultados[1]
    diferencia_loc_femalespcos.no.obese2[j] = vector_resultados[2]
    femalespcos.prueba.no.obese2[j] = vector_resultados[3]
    
    
    ### obesidad|grupo
    ### Vemos el sexo en las personas obesas
    
    Obese_var_i_sexo <- var_i[obesidad == "Obese" & grupo != "PCOS"]
    Obese_grupo_sexo <-  grupo[obesidad == "Obese" & grupo != "PCOS"]
    vector_resultados <-
      fn.paired(
        x = Obese_var_i_sexo,
        g = Obese_grupo_sexo,
        termino1 = "Female",
        termino2 = "Male"
      )
    sexo.pvalue.obese2[j] = vector_resultados[1]
    diferencia_loc_sexo.obese2[j] = vector_resultados[2]
    sexo.prueba.obese2[j] = vector_resultados[3]
    ###### Vemos hombres vs mujeres PCOS obesas
    
    Obese_var_i_malespcos <-
      var_i[obesidad == "Obese" & grupo != "Female"]
    Obese_grupo_malespcos <-
      grupo[obesidad == "Obese" & grupo != "Female"]
    vector_resultados <-
      fn.paired(
        x = Obese_var_i_malespcos,
        g = Obese_grupo_malespcos,
        termino1 = "Male",
        termino2 = "PCOS"
      )
    malespcos.pvalue.obese2[j] = vector_resultados[1]
    diferencia_loc_malespcos.obese2[j] = vector_resultados[2]
    malespcos.prueba.obese2[j] = vector_resultados[3]
    
    ### obesidad en mujeres vs mujeres PCOS
    Obese_var_i_femalespcos <-
      var_i[obesidad == "Obese" & grupo != "Male"]
    Obese_grupo_femalespcos <-
      grupo[obesidad == "Obese" & grupo != "Male"]
    vector_resultados <-
      fn.paired(
        x = Obese_var_i_femalespcos,
        g = Obese_grupo_femalespcos,
        termino1 = "Female",
        termino2 = "PCOS"
      )
    femalespcos.pvalue.obese2[j] = vector_resultados[1]
    diferencia_loc_femalespcos.obese2[j] = vector_resultados[2]
    femalespcos.prueba.obese2[j] = vector_resultados[3]
    
    
    ### AHORA VEREMOS EL GRUPO|OBESIDAD
    
    sexo_var_i <- var_i[grupo != "PCOS"]
    sexo.obesidad <- obesidad[grupo != "PCOS"]
    
    fems_var_i <- var_i[grupo != "Male"]
    fems.obesidad <- obesidad[grupo != "Male"]
    
    pcosmales_var_i <- var_i[grupo != "Female"]
    pcosmales.obesidad <- obesidad[grupo != "Female"]
    
    vector_resultados <-
      fn.paired(
        x = sexo_var_i,
        g = sexo.obesidad,
        termino1 = "No Obese",
        termino2 = "Obese"
      )
    sexo.p.value.obesidad2[j] <- vector_resultados[1]
    sexo.diferencia.obesidad2[j] <- vector_resultados[2]
    sexo.prueba.obesidad2[j] <- vector_resultados[3]
    
    vector_resultados <-
      fn.paired(
        x = fems_var_i,
        g = fems.obesidad,
        termino1 = "No Obese",
        termino2 = "Obese"
      )
    fems.p.value.obesidad2[j] <- vector_resultados[1]
    fems.diferencia.obesidad2[j] <- vector_resultados[2]
    fems.prueba.obesidad2[j] <- vector_resultados[3]
    
    
    vector_resultados <-
      fn.paired(
        x = pcosmales_var_i,
        g = pcosmales.obesidad,
        termino1 = "No Obese",
        termino2 = "Obese"
      )
    pcosmales.p.value.obesidad2[j] <- vector_resultados[1]
    pcosmales.diferencia.obesidad2[j] <- vector_resultados[2]
    pcosmales.prueba.obesidad2[j] <- vector_resultados[3]
    
    
    
  }
  
  
  ### CREAMOS LAS MATRICES DE RESULTADOS
  
  ####### MODELO 1 ########
  #################   OBESIDAD
  
  glm_modelo_1 <- cbind(
    obesidad = s.obbesidad.mod1,
    grupo = s.grupo.mod1,
    interaccion = s.interaccion.mod1,
    tipo_prueba = s.prueba.mod1
  )
  glm_modelo_2 <- cbind(
    obesidad = s.obbesidad.mod2,
    grupo = s.grupo.mod2,
    interaccion = s.interaccion.mod2,
    tipo_prueba = s.prueba.mod2
  )
  
  OBESIDAD_MODELO_1 <- cbind(p.value = obesidad.pvalue ,
                             diferencia = diferencia_loc_obesidad,
                             tipo.prueba = obesidad.prueba)
  
  ################# GRUPO
  
  GRUPO_MODELO_1 <- cbind(
    mu_Female_PCOS = diferencia_loc_femalespcos,
    mu_SEXO = diferencia_loc_sexo,
    mu_Males_PCOS = diferencia_loc_malespcos,
    p_value_Females_PCOS = femalespcos.pvalue,
    p_value_SEXO = sexo.pvalue,
    p_value_Males_PCOS = malespcos.pvalue,
    sexo_prueba = sexo.prueba,
    femalespcos_prueba = femalespcos.prueba,
    malespcos_prueba = malespcos.prueba
  )
  
  
  ################# NO OBESIDAD,OBESIDAD|grupo
  
  GRUPO_MODELO_NO_OBESE_1 <-
    cbind(
      mu_Female_PCOS_No_Obese = diferencia_loc_femalespcos.no.obese,
      mu_SEXO_No_Obese = diferencia_loc_sexo.no.obese,
      mu_Males_PCOS_No_Obese = diferencia_loc_malespcos.no.obese,
      p_value_Females_PCOS_No_Obese = femalespcos.pvalue.no.obese,
      p_value_SEXO_No_Obese = sexo.pvalue.no.obese,
      p_value_Males_PCOS_No_Obese = malespcos.pvalue.no.obese,
      sexo_prueba_No_Obese = sexo.prueba.no.obese,
      femalespcos_prueba_No_Obese = femalespcos.prueba.no.obese,
      malespcos_prueba_No_Obese = malespcos.prueba.no.obese
    )
  GRUPO_MODELO_OBESE_1 <-
    cbind(
      mu_Female_PCOS_Obese = diferencia_loc_femalespcos.obese,
      mu_SEXO_Obese = diferencia_loc_sexo.obese,
      mu_Males_PCOS_Obese = diferencia_loc_malespcos.obese,
      p_value_Females_PCOS_Obese = femalespcos.pvalue.obese,
      p_value_SEXO_Obese = sexo.pvalue.obese,
      p_value_Males_PCOS_Obese = malespcos.pvalue.obese,
      sexo_prueba_Obese = sexo.prueba.obese,
      femalespcos_prueba_Obese = femalespcos.prueba.obese,
      malespcos_prueba_Obese = malespcos.prueba.obese
    )
  
  ################## grupo| OBESIDAD
  
  GRUPO_OBESIDAD_GRUPO_1 <-
    cbind(
      mu_sexo_obesidad = sexo.diferencia.obesidad,
      mu_females_obesidad = fems.diferencia.obesidad,
      mu_pcosmales_obesidad = pcosmales.diferencia.obesidad,
      pvalue_sexo_obesidad = sexo.p.value.obesidad,
      pvalue_females_obesidad = fems.p.value.obesidad,
      pvalue_pcosmales_obesidad = pcosmales.p.value.obesidad,
      prueba_sexo_obesidad = sexo.prueba.obesidad,
      prueba_fems_obesidad = fems.prueba.obesidad,
      prueba_pcosmales_obesidad = pcosmales.prueba.obesidad
    )
  
  
  ####### MODELO 2 ########
  #################   OBESIDAD
  OBESIDAD_MODELO_2 <- cbind(p.value = obesidad.pvalue2 ,
                             diferencia = diferencia_loc_obesidad2,
                             tipo.prueba = obesidad.prueba2)
  
  ################# GRUPO
  
  GRUPO_MODELO_2 <-
    cbind(
      mu_Female_PCOS = diferencia_loc_femalespcos2,
      mu_SEXO = diferencia_loc_sexo2,
      mu_Males_PCOS = diferencia_loc_malespcos2,
      p_value_Females_PCOS = femalespcos.pvalue2,
      p_value_SEXO = sexo.pvalue2,
      p_value_Males_PCOS = malespcos.pvalue2,
      sexo_prueba = sexo.prueba2,
      femalespcos_prueba = femalespcos.prueba2,
      malespcos_prueba = malespcos.prueba2
    )
  
  
  ################# NO OBESIDAD,OBESIDAD|grupo
  
  GRUPO_MODELO_NO_OBESE_2 <-
    cbind(
      mu_Female_PCOS_No_Obese = diferencia_loc_femalespcos.no.obese2,
      mu_SEXO_No_Obese = diferencia_loc_sexo.no.obese2,
      mu_Males_PCOS_No_Obese = diferencia_loc_malespcos.no.obese2,
      p_value_Females_PCOS_No_Obese = femalespcos.pvalue.no.obese2,
      p_value_SEXO_No_Obese = sexo.pvalue.no.obese2,
      p_value_Males_PCOS_No_Obese = malespcos.pvalue.no.obese2,
      sexo_prueba_No_Obese = sexo.prueba.no.obese2,
      femalespcos_prueba_No_Obese = femalespcos.prueba.no.obese2,
      malespcos_prueba_No_Obese = malespcos.prueba.no.obese2
    )
  GRUPO_MODELO_OBESE_2 <-
    cbind(
      mu_Female_PCOS_Obese = diferencia_loc_femalespcos.obese2,
      mu_SEXO_Obese = diferencia_loc_sexo.obese2,
      mu_Males_PCOS_Obese = diferencia_loc_malespcos.obese2,
      p_value_Females_PCOS_Obese = femalespcos.pvalue.obese2,
      p_value_SEXO_Obese = sexo.pvalue.obese2,
      p_value_Males_PCOS_Obese = malespcos.pvalue.obese2,
      sexo_prueba_Obese = sexo.prueba.obese2,
      femalespcos_prueba_Obese = femalespcos.prueba.obese2,
      malespcos_prueba_Obese = malespcos.prueba.obese2
    )
  
  ################## grupo| OBESIDAD
  
  GRUPO_OBESIDAD_GRUPO_1 <-
    cbind(
      mu_sexo_obesidad = sexo.diferencia.obesidad2,
      mu_females_obesidad = fems.diferencia.obesidad2,
      mu_pcosmales_obesidad = pcosmales.diferencia.obesidad2,
      pvalue_sexo_obesidad = sexo.p.value.obesidad2,
      pvalue_females_obesidad = fems.p.value.obesidad2,
      pvalue_pcosmales_obesidad = pcosmales.p.value.obesidad2,
      prueba_sexo_obesidad = sexo.prueba.obesidad2,
      prueba_fems_obesidad = fems.prueba.obesidad2,
      prueba_pcosmales_obesidad = pcosmales.prueba.obesidad2
    )
  
  
  
  ###### Vemos hombres vs mujeres PCOS no obesas
  
  return(
    list(
      glm_modelo_1 = glm_modelo_1,
      OBESIDAD_MODELO_1 = OBESIDAD_MODELO_1,
      GRUPO_MODELO_1 = GRUPO_MODELO_1,
      GRUPO_MODELO_NO_OBESE_1 = GRUPO_MODELO_NO_OBESE_1,
      GRUPO_MODELO_OBESE_1 = GRUPO_MODELO_OBESE_1,
      GRUPO_OBESIDAD_GRUPO_1 = GRUPO_OBESIDAD_GRUPO_1
    )
  )
  
  
  
}
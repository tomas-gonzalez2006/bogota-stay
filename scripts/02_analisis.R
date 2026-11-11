# Ejecutar desde la raíz de Parcial_3.Rproj.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_ALL", ""))
source("scripts/00_configuracion.R", encoding="UTF-8", local=TRUE)


# ============================================================================
# 1. Matrices y mapas iniciales
# ============================================================================
# Tres conceptualizaciones sobre la misma muestra, sin usar y para elegir parámetros.
library(spdep)
library(patchwork)
b <- cargar_analisis()
bm <- st_transform(b, crs_metrico)
ids <- b$upz_codigo
centros <- suppressWarnings(st_point_on_surface(bm))
xy <- st_coordinates(centros)
distancias <- as.matrix(dist(xy))
diag(distancias) <- Inf
rownames(distancias) <- colnames(distancias) <- ids

# W1: contigüidad reina. Solo las islas se enlazan a su vecino más cercano,
# de forma recíproca. Se conserva la muestra y se documenta cada enlace añadido.
nb_reina <- poly2nb(bm, queen = TRUE, row.names = ids, snap = .5)
a1 <- nb2mat(nb_reina, style = "B", zero.policy = TRUE)
islas <- which(rowSums(a1) == 0)
puentes <- data.frame(origen = character(), destino = character(), distancia_m = numeric())
for (i in islas) {
  j <- which.min(distancias[i, ])
  a1[i, j] <- a1[j, i] <- 1
  puentes <- rbind(puentes, data.frame(origen = ids[i], destino = ids[j], distancia_m = distancias[i, j]))
}
write_csv(puentes, "resultados/pesos/enlaces_islas_W1.csv")
# W2: unión simétrica de k=4 vecinos. Aproximación a interacción de proximidad.
nb_k <- knn2nb(knearneigh(xy, k = 4), row.names = ids, sym = TRUE)
a2 <- nb2mat(nb_k, style = "B", zero.policy = FALSE)
# W3: umbral mínimo común que deja al menos un vecino en todas las unidades.
# Los representantes interiores evitan centroides fuera de polígonos cóncavos.
umbral <- max(apply(distancias, 1, min)) * (1 + 1e-6)
a3 <- 1L * (distancias <= umbral)
diag(a3) <- 0
dimnames(a1) <- dimnames(a2) <- dimnames(a3) <- list(ids, ids)
adj <- list(W1 = a1, W2 = a2, W3 = a3)
listw <- lapply(adj, function(a) mat2listw(a, row.names = ids, style = "W", zero.policy = FALSE))
stopifnot(all(vapply(adj, function(a) all(rowSums(a) > 0) && all(diag(a) == 0) && all(a == t(a)), logical(1))))
comparacion <- do.call(rbind, lapply(combn(names(adj), 2, simplify = FALSE), function(nm) {
  data.frame(matriz_a = nm[1], matriz_b = nm[2], enlaces_distintos = sum(adj[[nm[1]]] != adj[[nm[2]]]) / 2,
             matrices_identicas = identical(adj[[nm[1]]], adj[[nm[2]]]))
}))
stopifnot(!any(comparacion$matrices_identicas))
resumen <- bind_rows(lapply(names(listw), function(w) {
  lw <- listw[[w]]; suma <- rowSums(listw2mat(lw)); nc <- n.comp.nb(lw$neighbours)
  data.frame(matriz = w, regla = switch(w, W1 = "Reina + vecino recíproco solo para islas", W2 = "k=4; unión simétrica", W3 = "Umbral de distancia entre representantes interiores"),
    parametro = switch(w, W1 = "snap=0,5 m", W2 = "k=4", W3 = paste0(fmt(umbral / 1000), " km")),
    n = nrow(b), enlaces = sum(adj[[w]]) / 2, vecinos_min = min(card(lw$neighbours)),
    vecinos_media = mean(card(lw$neighbours)), vecinos_max = max(card(lw$neighbours)),
    islas_originales = if (w == "W1") length(islas) else 0L, islas_finales = sum(card(lw$neighbours) == 0),
    componentes = nc$nc, suma_filas_min = min(suma), suma_filas_max = max(suma),
    normalizacion = "W: filas suman 1", orden_verificado = identical(attr(lw$neighbours, "region.id"), ids))
}))
write_csv(resumen, "resultados/pesos/resumen_matrices.csv")
write_csv(comparacion, "resultados/pesos/comparacion_matrices.csv")
saveRDS(list(listw = listw, adyacencia = adj, coordenadas = xy, umbral_m = umbral, puentes = puentes,
             md5_base = unname(tools::md5sum("datos_limpios/base_analitica_upz.rds"))), "resultados/pesos/pesos_espaciales.rds")
for (w in names(listw)) {
  m <- listw2mat(listw[[w]])
  stopifnot(max(abs(rowSums(m) - 1)) < 1e-12)
  write_csv(data.frame(upz_codigo = ids, m, check.names = FALSE), paste0("resultados/pesos/", w, ".csv"))
}

# Clases de y: mismos datos, cinco clases con cuantiles e intervalos iguales.
q <- quantile(b$sentimiento_promedio, seq(0, 1, .2), names = FALSE)
eq <- seq(min(b$sentimiento_promedio), max(b$sentimiento_promedio), length.out = 6)
stopifnot(length(unique(q)) == 6)
etiquetas_clase <- function(cortes) paste0(fmt(head(cortes, -1), 2), "–", fmt(tail(cortes, -1), 2))
b$y_cuantiles <- cut(b$sentimiento_promedio, q, include.lowest = TRUE, labels = etiquetas_clase(q))
b$y_intervalos <- cut(b$sentimiento_promedio, eq, include.lowest = TRUE, labels = etiquetas_clase(eq))
palq <- setNames(hcl.colors(5, "YlOrRd", rev = TRUE), levels(b$y_cuantiles))
pale <- setNames(hcl.colors(5, "YlOrRd", rev = TRUE), levels(b$y_intervalos))
extremos <- b$upz_codigo[c(which.min(b$sentimiento_promedio), which.max(b$sentimiento_promedio))]
p1 <- mapa_categorico(b, "y_cuantiles", "Sentimiento · cuantiles", palq, "NRC ES (0–100)", "Cinco clases: aproximadamente igual número de UPZ", extremos)
p2 <- mapa_categorico(b, "y_intervalos", "Sentimiento · intervalos iguales", pale, "NRC ES (0–100)", "Cinco clases: igual amplitud numérica", extremos)
guardar_grafico(p1, "01_sentimiento_cuantiles")
guardar_grafico(p2, "02_sentimiento_intervalos")
guardar_grafico((p1 + labs(caption=NULL)) + (p2 + labs(caption=NULL)) + plot_layout(guides = "keep") +
                 plot_annotation(caption=fuente_mapa), "03_comparacion_clasificaciones", 16, 8)
guardar_grafico(mapa_numerico(b, "rating_estrellas", "Estrellas de los mismos anuncios", "Estrellas (1–5)", "Media por UPZ; mismos anuncios usados en sentimiento"), "04_estrellas")
guardar_grafico(mapa_numerico(b, "log_densidad", "Concentración de la oferta", "log(1 + anuncios/km²)", "Variable de control X2; anuncios completos del corte"), "05_densidad")
guardar_grafico(mapa_numerico(b, "desacople_percentiles", "Desacople de posiciones relativas", "Puntos percentiles", "Percentil de estrellas menos percentil de sentimiento", limites=c(-100,100), divergente = TRUE), "06_desacople")

# Comparación directa sobre posiciones relativas de los mismos anuncios y UPZ.
# Los índices originales se mantienen para modelos; no se fuerza equivalencia cardinal.
b$clase_percentil_rating <- clase_percentil(b$percentil_rating)
b$clase_percentil_sentimiento <- clase_percentil(b$percentil_sentimiento)
rotulos <- intersect(c("UPZ41", "UPZ99"), b$upz_codigo)
pr <- mapa_categorico(b,"clase_percentil_rating","Estrellas · posición relativa",colores_comparables,
  "Percentil en 82 UPZ\n(0–100)","Mismos anuncios pareados · media de estrellas por UPZ",rotulos)
ps <- mapa_categorico(b,"clase_percentil_sentimiento","Sentimiento · posición relativa",colores_comparables,
  "Percentil en 82 UPZ\n(0–100)","Mismos anuncios pareados · media de NRC ES por UPZ",rotulos)
guardar_grafico(pr,"17_estrellas_percentiles")
guardar_grafico(ps,"18_sentimiento_percentiles")
guardar_grafico((pr+labs(caption=NULL))+(ps+labs(caption=NULL))+
  plot_layout(guides="collect")+plot_annotation(
    title="Estrellas y sentimiento con la misma escala cartográfica",
    caption=paste0(fuente_mapa,"\nMismos cortes, colores, territorio y muestra. El percentil compara posición, no equivalencia entre estrellas y NRC.")),
  "19_comparacion_estrellas_sentimiento",16,8)
write_csv(st_drop_geometry(b) %>% select(upz_codigo,upz_nombre_clean,n_pareados,n_resenas_validas,
  rating_estrellas,sentimiento_promedio,percentil_rating,percentil_sentimiento,
  clase_percentil_rating,clase_percentil_sentimiento,desacople_percentiles),
  "resultados/tablas/comparacion_cartografica_pareada.csv")
write_csv(data.frame(desde=head(cortes_comparables,-1),hasta=tail(cortes_comparables,-1),
  etiqueta=etiquetas_comparables,color=unname(colores_comparables),
  cierre="Derecha cerrada; primer intervalo incluye 0",unidad="Percentil por UPZ en muestra fija"),
  "resultados/tablas/escala_comparable.csv")
write_csv(data.frame(upz_codigo = ids, upz = b$upz_nombre_clean, sentimiento = b$sentimiento_promedio,
                      cuantiles = b$y_cuantiles, intervalos_iguales = b$y_intervalos), "resultados/tablas/clasificaciones_y.csv")
write_csv(data.frame(metodo = rep(c("Cuantiles", "Intervalos iguales"), each = 5),
                     desde = c(head(q, -1), head(eq, -1)), hasta = c(tail(q, -1), tail(eq, -1)),
                     n = c(as.integer(table(b$y_cuantiles)), as.integer(table(b$y_intervalos)))), "resultados/tablas/cortes_y.csv")

# Redes: rótulos con códigos y claves disponibles en las tablas/atlas.
redes <- lapply(names(adj), function(w) {
  ij <- which(upper.tri(adj[[w]]) & adj[[w]] > 0, arr.ind = TRUE)
  segmentos <- data.frame(x = xy[ij[,1],1], y = xy[ij[,1],2], xe = xy[ij[,2],1], ye = xy[ij[,2],2])
  ggplot() + geom_sf(data = bm, fill = "#edf2f5", color = "#b8c7d0", linewidth = .2) +
    geom_segment(data = segmentos, aes(x=x, y=y, xend=xe, yend=ye), color = "#228c8d", alpha = .4, linewidth = .4) +
    geom_sf(data = centros, color = "#173850", size = 1) +
    labs(title = paste0(w, " · ", resumen$parametro[resumen$matriz == w]),
         subtitle = paste0(resumen$regla[resumen$matriz == w], "\n", sum(adj[[w]]) / 2, " enlaces; sin unidades aisladas"), caption = fuente_mapa) + tema_mapa()
})
guardar_grafico(wrap_plots(redes, nrow = 1), "07_tres_vecindades", 17, 7)
cat("Matrices y mapas iniciales terminados.\n"); print(resumen)


# ============================================================================
# 2. Autocorrelación global y local
# ============================================================================
# Inferencia global por permutaciones; LISA condicional con W1 principal.
library(spdep)
library(patchwork)
b <- cargar_analisis(); pesos <- cargar_pesos(b); n <- nrow(b)
# Las permutaciones globales reubican valores entre UPZ, manteniendo W fija.
set.seed(semilla)
perms <- replicate(permutaciones, sample.int(n))
indices <- function(x, W) {
  z <- x - mean(x); s0 <- sum(W)
  c(Moran = (length(x) / s0) * sum(z * as.vector(W %*% z)) / sum(z^2),
    Geary = ((length(x) - 1) / (2 * s0)) * sum(W * outer(x, x, "-")^2) / sum(z^2),
    Getis_Ord_G = sum(x * as.vector(W %*% x)) / (sum(x)^2 - sum(x^2)))
}
globales <- list(); contador <- 0L
for (variable in c("sentimiento_promedio", "rating_estrellas")) {
  x <- b[[variable]]
  for (w in c("W1", "W2")) {
    W <- listw2mat(pesos$listw[[w]])
    obs <- indices(x, W)
    simulados <- vapply(seq_len(permutaciones), function(i) indices(x[perms[, i]], W), numeric(3))
    esperados <- c(-1 / (n - 1), 1, sum(W) / (n * (n - 1)))
    # Dos colas definidas como distancia absoluta al valor esperado de H0.
    p <- vapply(1:3, function(j) (1 + sum(abs(simulados[j, ] - esperados[j]) >= abs(obs[j] - esperados[j]))) /
                  (permutaciones + 1), numeric(1))
    z <- (obs - rowMeans(simulados)) / apply(simulados, 1, sd)
    direccion <- c(if (obs[1] > esperados[1]) "similitud" else "disimilitud",
                   if (obs[2] < 1) "similitud" else "disimilitud",
                   if (obs[3] > esperados[3]) "concentracion de valores altos" else "concentracion de valores bajos")
    contador <- contador + 1L
    globales[[contador]] <- data.frame(variable = variable, matriz = w, indice = names(obs), estadistico = unname(obs),
      esperado_H0 = esperados, z_permutaciones = unname(z), p_valor = p,
      alternativa = "Bilateral: indice distinto de su referencia bajo aleatoriedad espacial",
      metodo = paste0(permutaciones, " permutaciones globales; distancia a E0; +1"), alfa = alfa,
      conclusion = ifelse(p < alfa, paste0("Rechazar H0: ", direccion), "No rechazar H0 de aleatoriedad espacial"))
  }
}
globales <- bind_rows(globales)
stopifnot(nrow(globales) == 12, all(is.finite(globales$p_valor)))
write_csv(globales, "resultados/tablas/autocorrelacion_global_12.csv")

lw <- pesos$listw$W1
y <- setNames(b$sentimiento_promedio, b$upz_codigo)
z <- as.numeric(scale(y)); lagz <- lag.listw(lw, z, zero.policy = FALSE)
# Guardar resultados originales permite auditar las columnas exactas de spdep.
moran <- localmoran_perm(y, lw, nsim = permutaciones, alternative = "two.sided", iseed = semilla,
                         no_repeat_in_row = TRUE, zero.policy = FALSE)
geary <- localC_perm(y, lw, nsim = permutaciones, alternative = "two.sided", iseed = semilla,
                     no_repeat_in_row = TRUE, zero.policy = FALSE)
# Gi (sin auto-vecino): usa exactamente W1, no una matriz W1+I diferente.
getis <- localG_perm(y, lw, nsim = permutaciones, alternative = "two.sided", iseed = semilla,
                     no_repeat_in_row = TRUE, zero.policy = FALSE)
saveRDS(list(moran = moran, geary = geary, getis = getis), "resultados/modelos/lisa_original.rds")
extraer_p <- function(tabla) {
  j <- which(grepl("Sim$", colnames(tabla)) & !grepl("folded", colnames(tabla)))
  if (length(j) != 1) stop("No se pudo identificar p bilateral de permutaciones: ", paste(colnames(tabla), collapse=", "))
  # Evita p=0 por discretización numérica del rango de simulaciones.
  pmax(1 / (permutaciones + 1), pmin(1, tabla[, j]))
}
p_m <- extraer_p(moran)
p_c <- extraer_p(attr(geary, "pseudo-p"))
p_g <- extraer_p(attr(getis, "internals"))
local <- st_drop_geometry(b) %>% select(upz_codigo, upz_nombre_clean, sentimiento_promedio) %>%
  mutate(z_y = z, lag_z_y = as.numeric(lagz), moran_Ii = moran[, "Ii"], p_moran = p_m,
         q_moran_BH = p.adjust(p_m, method = "BH"), geary_Ci = as.numeric(geary), p_geary = p_c,
         getis_Gi_z = as.numeric(getis), p_getis = p_g, q_getis_BH = p.adjust(p_g, method = "BH"),
         cuadrante = case_when(z_y >= 0 & lag_z_y >= 0 ~ "Alto–Alto", z_y < 0 & lag_z_y < 0 ~ "Bajo–Bajo",
                               z_y >= 0 & lag_z_y < 0 ~ "Alto–Bajo", TRUE ~ "Bajo–Alto"),
         moran_categoria = if_else(p_moran < alfa, cuadrante, "No significativo"),
         getis_categoria = case_when(p_getis >= alfa ~ "No significativo", getis_Gi_z > 0 ~ "Concentracion alta",
                                     TRUE ~ "Concentracion baja"))
write_csv(local, "resultados/tablas/autocorrelacion_local_W1.csv")
bl <- b %>% left_join(local %>% select(-upz_nombre_clean, -sentimiento_promedio), by = "upz_codigo", relationship = "one-to-one")
col_m <- c("Alto–Alto" = "#bb3846", "Bajo–Bajo" = "#286aa6", "Alto–Bajo" = "#f6a854", "Bajo–Alto" = "#7dbbda", "No significativo" = "#f3f0e9")
col_g <- c("Concentracion alta" = "#b52e43", "Concentracion baja" = "#2474a3", "No significativo" = "#f3f0e9")
ejemplos <- local$upz_codigo[order(local$p_moran)][1:3]
p1 <- mapa_categorico(bl, "moran_categoria", "Moran local · W1", col_m, "Categoría LISA",
                      "p bilateral < 0,05 · 9.999 permutaciones condicionales · sin ajuste múltiple", ejemplos)
p2 <- mapa_numerico(bl, "geary_Ci", "Geary local · W1", "C local (y tipificada)",
                    "Valores altos: mayor disimilitud con vecinos; bajos: mayor semejanza",
                    local$upz_codigo[c(which.max(local$geary_Ci), which.min(local$geary_Ci))])
p3 <- mapa_categorico(bl, "getis_categoria", "Getis–Ord Gi local · W1", col_g, "Concentración local",
                      "Gi sin auto-vecino · p bilateral < 0,05 · 9.999 permutaciones", local$upz_codigo[order(local$p_getis)][1:3])
guardar_grafico(p1, "08_moran_local"); guardar_grafico(p2, "09_geary_local"); guardar_grafico(p3, "10_getis_local")
guardar_grafico((p1 + labs(caption=NULL)) + (p2 + labs(caption=NULL)) + (p3 + labs(caption=NULL)) +
                 plot_annotation(caption=fuente_mapa), "11_comparacion_local", 18, 8)
png_moran <- ggplot(local, aes(z_y, lag_z_y)) + geom_hline(yintercept=0, color="#acb6bf") +
  geom_vline(xintercept=0, color="#acb6bf") + geom_point(aes(color=moran_categoria), size=2) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, color="#24364b", linewidth=.6) +
  scale_color_manual(values=col_m) + labs(title="Diagrama de Moran · sentimiento", x="Sentimiento tipificado (z)",
    y="Rezago espacial de z con W1", color="LISA", caption="Las categorías requieren p local < 0,05. El cuadrante por sí solo no implica significancia.") + theme_minimal()
guardar_grafico(png_moran, "12_diagrama_moran", 10, 6)
saveRDS(bl, "resultados/modelos/base_lisa.rds")
cat("12 pruebas globales y tres análisis locales terminados.\n")


# ============================================================================
# 3. Regresión: 16 ajustes y 18 diagnósticos
# ============================================================================
# Un MCO y 15 ajustes espaciales. Diagnóstico sobre innovaciones del modelo.
library(spdep)
library(spatialreg)
b <- cargar_analisis(); pesos <- cargar_pesos(b)
datos <- st_drop_geometry(b)
rownames(datos) <- b$upz_codigo
f <- sentimiento_promedio ~ rating_estrellas + log_densidad
y <- datos$sentimiento_promedio
modelos <- list(MCO = lm(f, data = datos, na.action = na.fail))
avisos <- data.frame(modelo=character(), aviso=character())
registrar <- function(expr, nombre) {
  withCallingHandlers(expr, warning=function(w) {
    avisos <<- rbind(avisos, data.frame(modelo=nombre, aviso=conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
}
for (w in names(pesos$listw)) {
  lw <- pesos$listw[[w]]
  dd <- datos
  dd$lag_rating_estrellas <- lag.listw(lw, dd$rating_estrellas)
  dd$lag_log_densidad <- lag.listw(lw, dd$log_densidad)
  modelos[[paste0("SLM_",w)]] <- registrar(lagsarlm(f, datos, lw, method="eigen", zero.policy=FALSE, na.action=na.fail), paste0("SLM_",w))
  modelos[[paste0("SEM_",w)]] <- registrar(errorsarlm(f, datos, lw, method="eigen", zero.policy=FALSE, na.action=na.fail), paste0("SEM_",w))
  modelos[[paste0("SLX_",w)]] <- lm(sentimiento_promedio ~ rating_estrellas + log_densidad + lag_rating_estrellas + lag_log_densidad,
                                      data=dd, na.action=na.fail)
  modelos[[paste0("SDM_",w)]] <- registrar(lagsarlm(f, datos, lw, Durbin=TRUE, method="eigen", zero.policy=FALSE, na.action=na.fail), paste0("SDM_",w))
  # La misma matriz entra en Wy y Wu.
  modelos[[paste0("SARAR_",w)]] <- registrar(sacsarlm(f, datos, lw, listw2=lw, method="eigen", zero.policy=FALSE,
                                                      na.action=na.fail, control=list(npars=9L)), paste0("SARAR_",w))
  cat("Estimadas cinco familias con", w, "\n")
}
stopifnot(length(modelos) == 16)
saveRDS(ajuste_compacto(modelos), "resultados/modelos/ajustes_16.rds")
write_csv(avisos, "resultados/modelos/avisos_estimacion.csv")

coeficientes <- list(); diag <- list(); residuos <- list(); impactos <- list(); k <- 0L
for (w in names(pesos$listw)) {
  W <- listw2mat(pesos$listw[[w]])
  for (familia in c("MCO", "SLM", "SEM", "SLX", "SDM", "SARAR")) {
    clave <- if (familia=="MCO") "MCO" else paste0(familia,"_",w)
    fit <- modelos[[clave]]; k <- k+1L
    if (inherits(fit,"lm")) {
      cf <- summary(fit)$coefficients
      cf_tab <- data.frame(termino=rownames(cf), estimacion=cf[,1], error_estandar=cf[,2], estadistico=cf[,3], p_valor=cf[,4], prueba="t de Student")
    } else {
      s <- summary(fit); cf <- s$Coef
      cf_tab <- data.frame(termino=rownames(cf), estimacion=cf[,1], error_estandar=cf[,2], estadistico=cf[,3], p_valor=cf[,4], prueba="Wald z asintotico")
      for (par in intersect(c("rho","lambda"), names(fit))) {
        se <- s[[paste0(par,".se")]]
        if (is.null(se)) se <- sqrt(diag(vcov(fit)))[par]
        z <- fit[[par]]/se
        cf_tab <- rbind(cf_tab, data.frame(termino=par, estimacion=fit[[par]], error_estandar=se, estadistico=z,
                                         p_valor=2*pnorm(-abs(z)), prueba="Wald z asintotico"))
      }
    }
    cf_tab$termino <- sub("^lag\\.", "lag_", cf_tab$termino)
    cf_tab$modelo <- familia; cf_tab$matriz <- w
    coeficientes[[k]] <- cf_tab
    # Construcción explícita de errores: SEM/SARAR se filtran por (I-lambda W).
    co <- coef(fit)
    if (inherits(fit,"Sarlm")) co <- fit$coefficients
    xb <- co["(Intercept)"] + co["rating_estrellas"] * datos$rating_estrellas + co["log_densidad"] * datos$log_densidad
    if (familia %in% c("SLX","SDM")) {
      nombres <- names(co)
      t1 <- nombres[grepl("^lag[._]rating_estrellas$", nombres)]
      t2 <- nombres[grepl("^lag[._]log_densidad$", nombres)]
      stopifnot(length(t1)==1, length(t2)==1)
      xb <- xb + co[t1] * as.vector(W %*% datos$rating_estrellas) + co[t2] * as.vector(W %*% datos$log_densidad)
    }
    rho <- if (familia %in% c("SLM","SDM","SARAR")) fit$rho else 0
    lambda <- if (familia %in% c("SEM","SARAR")) fit$lambda else 0
    error_estructural <- y - rho * as.vector(W %*% y) - xb
    e <- as.vector(error_estructural - lambda * as.vector(W %*% error_estructural))
    stopifnot(length(e)==nrow(b), all(is.finite(e)))
    # I residual con W del ajuste. Permutación descriptiva de residuos ajustados:
    # no es un test de especificación con parámetros reestimados bajo cada réplica.
    set.seed(semilla)
    mt <- moran.mc(e, pesos$listw[[w]], nsim=permutaciones, alternative="two.sided", zero.policy=FALSE)
    diag[[k]] <- data.frame(modelo=familia, matriz=w, n=length(e),
      pseudo_R2=1-var(e)/var(y), varianza_residual=var(e), moran_residual=as.numeric(mt$statistic),
      p_moran=mt$p.value, AIC=AIC(fit), logLik=as.numeric(logLik(fit)),
      tipo_residuo=if (lambda!=0) "Innovacion: (I-lambda W)(y-rho Wy-Xb-WXtheta)" else "Innovacion: y-rho Wy-Xb-WXtheta",
      inferencia_moran=paste0(permutaciones," permutaciones bilaterales de residuos; diagnostico exploratorio"))
    residuos[[paste0(familia,"_",w)]] <- e
    # Impactos exactos: diagonal media y suma media de filas del multiplicador.
    mult <- solve(diag(nrow(b))-rho*W)
    for (variable in c("rating_estrellas","log_densidad")) {
      beta <- co[variable]; theta <- 0
      if (familia %in% c("SLX","SDM")) theta <- co[grep(paste0("^lag[._]",variable,"$"), names(co))]
      S <- mult %*% (diag(nrow(b))*beta + W*theta)
      impactos[[length(impactos)+1L]] <- data.frame(modelo=familia, matriz=w, variable=variable,
        directo=mean(diag(S)), indirecto=mean(rowSums(S))-mean(diag(S)), total=mean(rowSums(S)))
    }
  }
}
diagnostico <- bind_rows(diag)
ct <- bind_rows(coeficientes)
stopifnot(nrow(diagnostico)==18, all(is.finite(diagnostico$AIC)), all(is.finite(ct$p_valor)))
write_csv(diagnostico, "resultados/tablas/diagnostico_18.csv")
write_csv(ct, "resultados/tablas/coeficientes_completos.csv")
# Tabla de los 16 ajustes sin repetir MCO; la versión de 18 filas facilita diagnóstico por W.
write_csv(ct %>% filter(modelo!="MCO" | matriz=="W1"), "resultados/tablas/coeficientes_16_ajustes.csv")
write_csv(bind_rows(impactos), "resultados/tablas/impactos_exactos.csv")
saveRDS(residuos, "resultados/modelos/residuos.rds")
capture.output(lapply(modelos, summary), file="resultados/modelos/resumenes_16.txt")
cat("16 ajustes, 18 diagnósticos y coeficientes completos terminados.\n")


# ============================================================================
# 4. Robustez y mapas residuales
# ============================================================================
library(patchwork)
b <- cargar_analisis(); pesos <- cargar_pesos(b)
diagnostico <- read_csv("resultados/tablas/diagnostico_18.csv", show_col_types=FALSE)
coefs <- read_csv("resultados/tablas/coeficientes_completos.csv", show_col_types=FALSE)
residuos <- readRDS("resultados/modelos/residuos.rds")
ranking <- diagnostico %>% group_by(matriz) %>% arrange(AIC, .by_group=TRUE) %>%
  mutate(orden_AIC=row_number(), delta_AIC=AIC-min(AIC)) %>% ungroup()
write_csv(ranking, "resultados/tablas/ranking_modelos.csv")
ganadores <- ranking %>% filter(orden_AIC==1)
write_csv(ganadores, "resultados/tablas/menor_AIC_por_W.csv")
# Comparación de las dos X centrales con familia fija: conservar signos y p concretos.
robustez <- coefs %>% filter(termino %in% c("rating_estrellas","log_densidad")) %>%
  mutate(significativo=p_valor<alfa) %>% group_by(modelo, termino) %>% summarise(
    beta_W1=estimacion[matriz=="W1"], p_W1=p_valor[matriz=="W1"],
    beta_W2=estimacion[matriz=="W2"], p_W2=p_valor[matriz=="W2"],
    beta_W3=estimacion[matriz=="W3"], p_W3=p_valor[matriz=="W3"],
    signo_estable=n_distinct(sign(estimacion))==1,
    significancia_estable=n_distinct(significativo)==1,
    rango_magnitud=max(estimacion)-min(estimacion), .groups="drop") %>%
  mutate(clasificacion=if_else(signo_estable & significancia_estable,
                               "Signo y decisión p estables (revisar magnitudes)", "Sensible a la matriz"))
write_csv(robustez, "resultados/tablas/robustez_coeficientes.csv")

# Elegir con W1 por fundamento territorial; delta AIC <2 favorece parsimonia MCO.
orden_w1 <- ranking %>% filter(matriz=="W1") %>% arrange(AIC)
elegido <- if (orden_w1$delta_AIC[orden_w1$modelo=="MCO"] < 2) "MCO" else orden_w1$modelo[1]
# La misma familia espacial se representa con cada W incluso si se prefiere MCO.
familia_mapa <- if (elegido=="MCO") orden_w1$modelo[orden_w1$modelo!="MCO"][1] else elegido
seleccion <- data.frame(modelo=elegido, matriz="W1", familia_residuos=familia_mapa,
  criterio="W1 por continuidad territorial; AIC junto con coeficientes y Moran residual; MCO preferido si delta AIC <2")
write_csv(seleccion, "resultados/tablas/seleccion_modelo.csv")

claves <- c("MCO_W1", paste0(familia_mapa,"_",c("W1","W2","W3")))
maximo <- max(abs(unlist(residuos[claves])))
cortes <- seq(-maximo,maximo,length.out=7)
labels <- paste0(fmt(head(cortes,-1),2)," a ",fmt(tail(cortes,-1),2))
colores <- setNames(c("#164c78","#4c91bb","#bedde9","#f3d8c8","#d97762","#92283d"), labels)
mapas <- lapply(claves,function(clave) {
  bb <- b; bb$residuo <- residuos[[clave]]
  bb$clase_residuo <- cut(bb$residuo,cortes,include.lowest=TRUE,labels=labels)
  met <- diagnostico[paste0(diagnostico$modelo,"_",diagnostico$matriz)==clave,]
  titulo <- if (clave=="MCO_W1") "MCO · ajuste único" else sub("_"," · ",clave)
  p <- mapa_categorico(bb,"clase_residuo",titulo,colores,"Innovación\n(puntos de NRC)",
    paste0("Cortes comunes · I residual = ",fmt(met$moran_residual),"; p = ",fp(met$p_moran)))
  guardar_grafico(p,paste0("13_residuos_",clave))
  p
})
mapas_panel <- lapply(mapas,function(p) p + labs(caption=NULL) + theme(plot.title=element_text(size=14),plot.subtitle=element_text(size=9)))
guardar_grafico(wrap_plots(mapas_panel,ncol=2,guides="collect") + plot_annotation(
  title="MCO y la misma familia espacial con tres matrices",caption=fuente_mapa) & theme(legend.position="right"),
  "14_cuatro_residuos", 14, 13)
write_csv(data.frame(desde=head(cortes,-1),hasta=tail(cortes,-1)), "resultados/tablas/cortes_residuos.csv")
res_df <- data.frame(upz_codigo=b$upz_codigo,upz=b$upz_nombre_clean,as.data.frame(residuos),check.names=FALSE)
write_csv(res_df,"resultados/tablas/residuos_18.csv")

# Magnitud y patrón de cambios espaciales, con nombres y valores concretos.
comparacion_res <- bind_rows(lapply(claves[-1],function(clave) {
  data.frame(upz_codigo=b$upz_codigo,upz=b$upz_nombre_clean,modelo_matriz=clave,
    residuo_MCO=residuos$MCO_W1,residuo_espacial=residuos[[clave]],
    cambio=residuos[[clave]]-residuos$MCO_W1,
    conserva_signo=sign(residuos[[clave]])==sign(residuos$MCO_W1))
}))
write_csv(comparacion_res,"resultados/tablas/cambios_residuos.csv")
graf <- coefs %>% filter(termino %in% c("rating_estrellas","log_densidad")) %>%
  mutate(termino=factor(termino,levels=c("rating_estrellas","log_densidad"),labels=c("Estrellas","log(1 + densidad)")),
         modelo=factor(modelo,levels=rev(c("MCO","SLM","SEM","SLX","SDM","SARAR"))))
p <- ggplot(graf,aes(estimacion,modelo,color=matriz)) + geom_vline(xintercept=0,color="#c2c9ce") +
  geom_errorbar(aes(xmin=estimacion-1.96*error_estandar,xmax=estimacion+1.96*error_estandar),
                 position=position_dodge(width=.5),width=.15,orientation="y") +
  geom_point(position=position_dodge(width=.5),size=2) + facet_wrap(~termino,scales="free_x") +
  scale_color_manual(values=c(W1="#165d85",W2="#ce6d3a",W3="#3c8b6d")) + theme_minimal(base_size=12) +
  labs(title="Coeficientes de las dos explicativas centrales",x="Coeficiente e intervalo Wald aproximado 95 %",y=NULL,
       caption="Intervalos clásicos/asintóticos. En modelos con Wy, los coeficientes requieren interpretación mediante impactos.\nMCO es el mismo ajuste en W1–W3; no se ha reestimado.")
guardar_grafico(p,"15_coeficientes_robustez",12,7)
p <- ggplot(ranking,aes(modelo,delta_AIC,fill=matriz)) + geom_col(position="dodge") +
  geom_hline(yintercept=2,linetype="dashed",color="#686e78") + theme_minimal(base_size=12) +
  scale_fill_manual(values=c(W1="#165d85",W2="#ce6d3a",W3="#3c8b6d")) +
  labs(title="Comparación de ajuste penalizado",x=NULL,y="Diferencia de AIC respecto del mínimo de cada W",caption="Misma respuesta y 82 UPZ en las 18 combinaciones.")
guardar_grafico(p,"16_comparacion_AIC",10,6)
cat("Tablas de robustez y cuatro mapas de residuos terminados. Selección:",elegido,"/ W1; contraste:",familia_mapa,"\n")


# ============================================================================
# 5. Sensibilidad e influencia
# ============================================================================
# Diagnósticos complementarios. No modifican la muestra de los 16 ajustes centrales.
library(lmtest)
library(sandwich)
b <- cargar_analisis(); fits <- readRDS("resultados/modelos/ajustes_16.rds")
m <- fits$MCO; e <- residuals(m)
hc3 <- coeftest(m,vcov.=vcovHC(m,type="HC3"))
write_csv(data.frame(termino=rownames(hc3),estimacion=hc3[,1],error_HC3=hc3[,2],t=hc3[,3],p_HC3=hc3[,4]),
          "resultados/tablas/MCO_HC3.csv")
bp <- bptest(m)
sh <- shapiro.test(e)
write_csv(data.frame(prueba=c("Breusch-Pagan","Shapiro-Wilk"),estadistico=c(unname(bp$statistic),unname(sh$statistic)),
                     p_valor=c(bp$p.value,sh$p.value),
                     lectura=c("H0: varianza constante","H0: residuos normales")),"resultados/tablas/diagnostico_adicional_MCO.csv")
influencia <- data.frame(upz_codigo=b$upz_codigo,upz=b$upz_nombre_clean,Cook=cooks.distance(m),
                         leverage=hatvalues(m),residuo_student=rstudent(m)) %>% arrange(desc(Cook))
write_csv(influencia,"resultados/tablas/influencia_MCO.csv")
loo <- bind_rows(lapply(seq_len(nrow(b)),function(i) {
  mm <- lm(sentimiento_promedio~rating_estrellas+log_densidad,data=st_drop_geometry(b)[-i,])
  cf <- coef(summary(mm))
  data.frame(upz_excluida=b$upz_codigo[i],nombre=b$upz_nombre_clean[i],beta_estrellas=cf["rating_estrellas",1],
             p_estrellas=cf["rating_estrellas",4],beta_densidad=cf["log_densidad",1],p_densidad=cf["log_densidad",4])
}))
write_csv(loo,"resultados/tablas/influencia_dejar_una_UPZ.csv")


# Paquete compacto: termina tras HC3 e influencia.
# Las sensibilidades por anuncio e idioma requieren el proyecto completo.

# Parámetros y funciones comunes; no ejecuta la preparación ni los modelos.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_ALL", ""))


# ============================================================================
# Parámetros del estudio y funciones de datos
# ============================================================================
# Parámetros compartidos: ejecutar desde Parcial_3.Rproj.
options(stringsAsFactors = FALSE)
fecha_inicio <- as.Date("2025-01-01")
fecha_fin <- as.Date("2026-06-22") # Último corte observado en listings.csv.
crs_metrico <- 3116L # MAGNA-SIRGAS / Colombia Bogota zone, metros.
min_alojamientos_pareados <- 3L
min_resenas_validas <- 10L
semilla <- 20261110L
alfa <- 0.05
permutaciones <- 9999L
for (carpeta in c("datos_limpios", "resultados", "resultados/calidad", "mapas")) {
  dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)
}

exigir_columnas <- function(x, columnas, archivo) {
  faltan <- setdiff(columnas, names(x))
  if (length(faltan)) stop(archivo, ": faltan columnas: ", paste(faltan, collapse = ", "))
}

# Se leen los IDs como texto: los IDs largos de Airbnb exceden la precisión de double.
leer_csv_texto <- function(ruta) {
  x <- readr::read_csv(ruta, col_types = readr::cols(.default = readr::col_character()),
                       na = c("", "NA", "N/A"), progress = FALSE)
  if (nrow(readr::problems(x))) stop("Problemas de lectura CSV en ", ruta)
  x
}

unicos_por_id <- function(x, campo, etiqueta) {
  originales <- nrow(x)
  x <- dplyr::distinct(x)
  if (anyNA(x[[campo]]) || anyDuplicated(x[[campo]])) {
    stop(etiqueta, ": IDs faltantes o duplicados con contenido diferente; resolver en la fuente.")
  }
  readr::write_csv(data.frame(archivo = etiqueta, filas_originales = originales,
                             duplicados_exactos_retirados = originales - nrow(x), filas = nrow(x)),
                   paste0("resultados/calidad/duplicados_", etiqueta, ".csv"))
  x
}

media_segura <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
mediana_segura <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}
ponderada_segura <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (any(ok)) sum(x[ok] * w[ok]) / sum(w[ok]) else NA_real_
}

# Las fórmulas solo usan columnas del modelo. Evitar serializar el entorno
# completo del análisis (mapas, permutaciones y objetos temporales) en cada ajuste.
ajuste_compacto <- function(x) {
  if (is.list(x)) {
    for (i in seq_along(x)) x[i] <- list(ajuste_compacto(x[[i]]))
  }
  if (!is.null(attr(x,".Environment",exact=TRUE))) attr(x,".Environment") <- baseenv()
  terminos <- attr(x,"terms",exact=TRUE)
  if (!is.null(terminos)) attr(x,"terms") <- ajuste_compacto(terminos)
  x
}


# ============================================================================
# Lectura de resultados, tablas y mapas PNG
# ============================================================================
# Utilidades compartidas por las etapas 04–09. Ejecutar desde la raíz del Rproj.
suppressPackageStartupMessages({library(sf); library(dplyr); library(readr); library(ggplot2)})
for (carpeta in c("resultados/tablas", "resultados/modelos", "resultados/pesos", "presentacion"))
  dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)

cargar_analisis <- function() {
  b <- readRDS("datos_limpios/base_analitica_upz.rds")
  b <- b[order(b$orden_unidad), ]
  class(b) <- c("sf", "data.frame")
  if (anyDuplicated(b$upz_codigo) || anyNA(st_drop_geometry(b)[c("sentimiento_promedio", "rating_estrellas", "log_densidad")]))
    stop("La muestra analítica tiene duplicados o faltantes.")
  rownames(b) <- b$upz_codigo
  b
}
cargar_pesos <- function(b) {
  p <- readRDS("resultados/pesos/pesos_espaciales.rds")
  for (w in p$listw) stopifnot(identical(attr(w$neighbours, "region.id"), b$upz_codigo))
  stopifnot(identical(p$md5_base, unname(tools::md5sum("datos_limpios/base_analitica_upz.rds"))))
  p
}
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
fp <- function(x) ifelse(is.na(x), "NA", ifelse(x < .001, "<0,001", gsub("\\.", ",", fmt(x, 3))))
# Contrato común de comparación cartográfica. Percentiles, no notas equivalentes.
cortes_comparables <- seq(0, 100, 20)
etiquetas_comparables <- c("0–20", "20–40", "40–60", "60–80", "80–100")
colores_comparables <- setNames(c("#f4efbd", "#c8dcb6", "#87b8be", "#4b89a7", "#174f70"), etiquetas_comparables)
clase_percentil <- function(x) cut(x, cortes_comparables, include.lowest=TRUE,
                                  right=TRUE, labels=etiquetas_comparables)
escape_html <- function(x) as.character(htmltools::htmlEscape(as.character(x)))
tabla_html <- function(x, digitos = 3) {
  z <- as.data.frame(x)
  for (j in seq_along(z)) {
    if (is.numeric(z[[j]])) z[[j]] <- ifelse(is.na(z[[j]]), "—", formatC(z[[j]], digits = digitos, format = "f"))
    else z[[j]] <- as.character(z[[j]])
  }
  paste0("<div class='tabla'><table><thead><tr>", paste0("<th>", escape_html(names(z)), "</th>", collapse = ""),
         "</tr></thead><tbody>", paste(apply(z, 1, function(r) paste0("<tr>", paste0("<td>", escape_html(r), "</td>", collapse = ""), "</tr>")), collapse = ""),
         "</tbody></table></div>")
}
fuente_mapa <- "Reseñas ES: 2025-01-01 a 2026-06-21 · Anuncios: junio de 2026\nFuentes: Inside Airbnb; UPZ (repositorio Duke/im90). Elaboración propia · EPSG:3116."
tema_mapa <- function() theme_void(base_size = 11) + theme(
  plot.title = element_text(face = "bold", color = "#172a3a", size = 16),
  plot.subtitle = element_text(color = "#485868", size = 10, margin = margin(b = 10)),
  plot.caption = element_text(hjust = 0, size = 8, color = "#53616f", margin = margin(t = 8)),
  legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 9),
  plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(14, 14, 14, 14))
guardar_grafico <- function(p, nombre, ancho = 10, alto = 8) {
  ggsave(paste0("mapas/", nombre, ".png"), p, width = ancho, height = alto, dpi = 180, bg = "white")
}
mapa_numerico <- function(b, variable, titulo, leyenda, subtitulo = "", etiquetas = NULL, limites = NULL, divergente = FALSE) {
  fondo <- st_transform(readRDS("datos_limpios/upz_agregado_analisis.rds"), crs_metrico)
  bm <- st_transform(b, crs_metrico)
  p <- ggplot() + geom_sf(data = fondo, fill = "#eceff1", color = "white", linewidth = .2) +
    geom_sf(data = bm, aes(fill = .data[[variable]]), color = "#ffffff", linewidth = .25)
  if (divergente) p <- p + scale_fill_gradient2(low = "#24629a", mid = "#faf7f0", high = "#bc3a41", midpoint = 0, limits = limites, name = leyenda)
  else p <- p + scale_fill_viridis_c(option = "C", direction = -1, limits = limites, name = leyenda)
  if (!is.null(etiquetas)) {
    et <- bm[bm$upz_codigo %in% etiquetas, ]
    if (nrow(et)) p <- p + geom_sf_label(data = suppressWarnings(st_point_on_surface(et)), aes(label = upz_nombre_clean),
                                        size = 2.5, linewidth = .15, alpha = .95, fill="white", color="#172a3a")
  }
  p + labs(title = titulo, subtitle = subtitulo, caption = paste0(fuente_mapa, "\nGris: UPZ fuera de la muestra; n = ", nrow(bm), ".")) + tema_mapa()
}
mapa_categorico <- function(b, variable, titulo, colores, leyenda, subtitulo = "", etiquetas = NULL) {
  fondo <- st_transform(readRDS("datos_limpios/upz_agregado_analisis.rds"), crs_metrico)
  bm <- st_transform(b, crs_metrico)
  p <- ggplot() + geom_sf(data = fondo, fill = "#e0e5e9", color = "white", linewidth = .2) +
    geom_sf(data = bm, aes(fill = .data[[variable]]), color = "white", linewidth = .25) +
    scale_fill_manual(values = colores, drop = FALSE, name = leyenda)
  if (!is.null(etiquetas)) {
    et <- bm[bm$upz_codigo %in% etiquetas, ]
    if (nrow(et)) p <- p + geom_sf_label(data = suppressWarnings(st_point_on_surface(et)), aes(label = upz_nombre_clean),
                                        size = 2.5, linewidth = .15, alpha = .95, fill="white", color="#172a3a")
  }
  p + labs(title = titulo, subtitle = subtitulo, caption = paste0(fuente_mapa, "\nGris de fondo: fuera de la muestra; n = ", nrow(bm), ".")) + tema_mapa()
}

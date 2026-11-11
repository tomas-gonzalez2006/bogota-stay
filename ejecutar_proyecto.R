if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_ALL", ""))
# Abrir Parcial_3.Rproj antes de ejecutar.
paquetes <- c("sf", "dplyr", "readr", "ggplot2", "spdep", "spatialreg", "patchwork", "lmtest", "sandwich")
instalar_dependencias <- function() {
  faltan <- setdiff(paquetes, rownames(installed.packages()))
  if (length(faltan)) install.packages(faltan, repos="https://cloud.r-project.org")
}
ejecutar_proyecto <- function() {
  if (!file.exists("datos_limpios/base_analitica_upz.rds")) stop("Abre Parcial_3.Rproj desde la carpeta descomprimida.")
  faltan <- setdiff(paquetes, rownames(installed.packages()))
  if (length(faltan)) stop("Ejecuta instalar_dependencias(). Faltan: ", paste(faltan, collapse=", "))
  source("scripts/02_analisis.R", encoding="UTF-8", local=new.env(parent=globalenv()))
  message("Listo: mapas/ y resultados/. El paquete reproduce el analisis desde las 82 UPZ preparadas.")
  invisible(TRUE)
}
demostracion_defensa <- function() {
  b <- readRDS("datos_limpios/base_analitica_upz.rds")
  m <- lm(sentimiento_promedio ~ rating_estrellas + log_densidad, data=sf::st_drop_geometry(b))
  print(coef(summary(m)))
  cat("Cambio por +0.1 estrellas:", 0.1 * coef(m)["rating_estrellas"], "puntos NRC\n")
  invisible(m)
}
if (sys.nframe() == 0L) ejecutar_proyecto()

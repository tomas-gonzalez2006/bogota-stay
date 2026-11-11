# Bogotá Stay: código para reutilizar

Autor: Tomás Alejandro González Forero. Análisis espacial, 2026-II.

## Ejecución

Guarda esta carpeta completa en tu computador. Con R y RStudio instalados, abre `Parcial_3.Rproj` y ejecuta:

```r
source("ejecutar_proyecto.R", encoding = "UTF-8")
instalar_dependencias()  # Primera vez; necesita Internet.
ejecutar_proyecto()
```

Genera mapas PNG en `mapas/` y tablas, pesos y modelos en `resultados/`.
Para un cálculo breve, después de instalar dependencias:

```r
mco <- demostracion_defensa()
```

Probado con R 4.4.3 en Windows. En Linux, sf puede necesitar bibliotecas de sistema. Las 9.999 permutaciones pueden tardar varios minutos. No requiere Python ni descargas de datos; después de instalar paquetes funciona sin Internet.

## Contenido mínimo

- `ejecutar_proyecto.R`: instalación, ejecución y demostración MCO.
- `scripts/00_configuracion.R`: parámetros y funciones comunes.
- `scripts/02_analisis.R`: tres matrices, 12 pruebas globales, tres indicadores locales, 16 modelos, diagnósticos, robustez entre W, mapas, HC3 e influencia.
- `datos_limpios/base_analitica_upz.rds`: base de 82 UPZ con geometría; entrada de los scripts.
- `datos_limpios/base_analitica_upz.csv`: la base sin geometría, para explorar en otras herramientas.
- `datos_limpios/upz_agregado_analisis.rds`: las 112 UPZ para dibujar el fondo de los mapas; los modelos usan solo las 82 de la base analítica.
- `docs/DICCIONARIO_DATOS.md`: definición de columnas.
- `Parcial_3.Rproj`: acceso desde RStudio.
- `docs/versiones_paquetes.csv`: versiones del entorno original como referencia; el instalador usa las disponibles en CRAN.

Edición compacta del proyecto original. El análisis conserva el código original hasta el diagnóstico de influencia. No incluye preparación de reseñas/NLP, sensibilidades por idioma/cobertura/borde ni generación de diapositivas, informe o atlas. Estas etapas requieren el proyecto completo. La demostración compacta reproduce MCO; la del proyecto completo consulta además W1 y Moran local.

## Datos e interpretación

Una fila es una UPZ: 82 UPZ, 9.625 anuncios pareados y 127.308 reseñas válidas en español. Reseñas: 01/01/2025–21/06/2026. Corte de anuncios: 21–22/06/2026.

Respuesta: sentimiento NRC en español (0–100). Explicativas: estrellas acumuladas (1–5) y log(1 + anuncios/km²). Las estrellas históricas y el texto reciente no miden exactamente el mismo periodo. La asociación no es causal y la inferencia de estrellas es sensible a influencia y errores robustos.

El RDS conserva geometría y orden de unidades. Para otro territorio, adapta columnas, rótulos y vecindades: varias etiquetas son específicas de Bogotá. El diccionario describe la base original; esta edición solo incorpora los formatos RDS y CSV.

## Fuentes

- Datos: Inside Airbnb, Bogotá, corte local junio de 2026: https://insideairbnb.com/ . No se conservó el enlace histórico exacto.
- Cartografía UPZ: im90, repositorio académico Duke: https://gitlab.oit.duke.edu/im90/enriched-spatial-analyses-of-pm2.5/-/raw/main/shapefiles/upz.geojson . No se afirma vigencia administrativa como UPL.
- Medición: Mohammad y Turney (2013), https://doi.org/10.1111/j.1467-8640.2012.00460.x .
- Contexto: Zervas et al. (2021), https://doi.org/10.1007/s11002-020-09546-4 ; Gurran y Phibbs (2017), https://doi.org/10.1080/01944363.2016.1249011 .
- Métodos: Anselin (1995), https://doi.org/10.1111/j.1538-4632.1995.tb00338.x ; Bivand y Wong (2018), https://doi.org/10.1007/s11749-018-0599-x .

Al reutilizar, reconoce al autor y las fuentes de datos y métodos. El paquete no añade una licencia a datos de terceros.

Verificación: esta edición se ejecutó desde las bases incluidas y generó las matrices, 12 pruebas globales, indicadores locales y 16 ajustes. Las versiones instaladas en otro equipo pueden producir pequeñas diferencias numéricas.

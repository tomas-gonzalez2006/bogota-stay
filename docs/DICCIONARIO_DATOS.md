# Diccionario de la base analítica

Archivo principal: `datos_limpios/base_analitica_upz.rds`. Una fila por UPZ incluida en la muestra. El CSV contiene las mismas columnas sin geometría. El RDS `upz_agregado_analisis.rds` incluye las 112 UPZ para el fondo cartográfico; los modelos usan solo las 82 de la base analítica. El `README.md` explica la pregunta, ejecución y alcance.

| Campo | Definición |
|---|---|
| `upz_codigo` | Identificador territorial estable de la capa original. Texto. |
| `upz_nombre_clean` | Nombre real de la UPZ; no OBJECTID. |
| `localidad_clean` | Localidad(es) de las partes de la UPZ en la fuente. No se escoge arbitrariamente una si existen varias. |
| `orden_unidad` | Orden alfabético por código para alinear datos y W; no es un ID sustantivo. |
| `area_km2` | Área de la UPZ disuelta, calculada en EPSG:3116 y convertida a km². |
| `n_alojamientos` | Todos los anuncios completos con coordenadas válidas y asignación única a esa UPZ. |
| `total_resenas` | Conteo de reseñas del periodo para esos anuncios, sin restricción de idioma. No es ocupación ni número de huéspedes. |
| `n_pareados` | Anuncios con estrellas válidas y al menos una reseña en español con cobertura NRC. |
| `n_resenas_validas` | Reseñas en español con cobertura NRC de los anuncios pareados. |
| `cobertura_pareados` | `n_pareados / n_alojamientos`. |
| `sentimiento_promedio` | **y**: promedio simple de los sentimientos medios por anuncio, usando español. Índice 0–100. |
| `rating_estrellas` | **X1**: promedio simple de estrellas de esos mismos anuncios, escala 1–5. |
| `rating_promedio` | Transformación visual: 20 × `rating_estrellas`. No iguala el constructo al sentimiento. |
| `densidad_alojamientos` | `n_alojamientos / area_km2`; anuncios por km². |
| `log_densidad` | **X2**: logaritmo natural de 1 + densidad. |
| `percentil_rating`, `percentil_sentimiento` | Rangos medios convertidos a percentiles en la muestra fija; empates reciben rango promedio. |
| `desacople_percentiles` | Percentil de rating menos percentil de sentimiento. Positivo: mejor posición relativa por estrellas. |
| `rating_ponderado_resenas` | Rating 0–100 ponderado por reseñas válidas en español de cada anuncio pareado. Auxiliar. |
| `sentimiento_ponderado_resenas` | Sentimiento ponderado por los mismos conteos. Auxiliar. |
| `sentimiento_es` | Media zonal de sentimiento en español; coincide con y. |
| `sentimiento_en` | Media de sentimiento en inglés entre los anuncios pareados que además tienen inglés válido. Puede usar un subconjunto. Auxiliar, no comparación plenamente pareada por idioma. |
| `sentimiento_bilingue_auxiliar` | Media del índice bilingüe de los anuncios pareados. No usar como y principal. |
| `n_listings_es`, `n_listings_en` | Número de anuncios pareados con sentimiento válido en cada idioma. |
| `precio_mediano_auxiliar` | Mediana de precios positivos finitos, sin excluir atípicos marcados. Moneda/unidad por confirmar con metadatos; no se usa como X. |
| `n_cerca_borde` | Anuncios de la oferta a ≤150 m del borde de su UPZ. |
| `motivo_muestra` | En la base final es `incluida`; en la base zonal completa explica exclusiones. |
| Geometría sf | Polígonos disueltos, EPSG:4326 en exportaciones. |


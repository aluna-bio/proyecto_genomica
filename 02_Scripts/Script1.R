#Instalar paquetes a utilizar

BiocManager::install("WGCNA", force = TRUE)
library(WGCNA)
install.packages("igraph")
library(igraph)

# Investigamos y esta función ignora los demás metadatos de la matriz, dejando solo los datos de expresión
#datos <- read.delim("01_RawData/GSE25066_series_matrix.txt", comment.char = "!", header = TRUE)

# Pasamos la primera columna, ID de los genes a los nombres de las filas (rownames)
#rownames(datos) <- datos[, 1]

# Eliminamos la columna de los ID
#datos <- datos [,-1]

# Calculamos la varianza de toda la matriz fila por fila 
# (¿que genes presentan la mayor diferencia entre muestras?) 
# (varianza=0 no hay diferencias entre las mjuestras)
#varianzas <- apply(datos,1,var)

# Ordenamos los genes de mayor varianza a menor varianza
# Nos quedamos con los primeros 5,000 - SOLO CON LOS NOMBRES
# Los guardamos en una lista -> genes_variables
#genes_variables <- names(sort(varianzas, decreasing = TRUE)[1:5000])

# De la matriz original seleccionamos las filas de los 5,000 genes mas variables
#datos_filtrados <- datos[genes_variables, ]

# Lo guardamos en un objeto RDS para evitar el procesamineto previo ya que la base de datos es muy pesada para Github, en el README se encuentra el link a la base de datos original.
#saveRDS(datos_filtrados, file = "01_RawData/datos.rds") # Este archivo se encuenra guardado en la carpeta 01_RawData

# Leemos el archivo RDS y lo guardamos en un objeto
datos_filtrados <- readRDS( file = "01_RawData/datos.rds" )

# Realizamos la transposición de la matriz, para cambiar los genes de fila a columna y las muestras de columna a fila
datos_t <- t(datos_filtrados)

# Busqueda del power optimo.
# Determinar el numero por el cual se potenciara toda la matriz
# Esto fortalece la correlación fuerte y debilita la correlación débil o negativa, haciendo que el valor de correlación sea más consistente con Red free scale. Las características son biológicamente más significativas. 
# (https://programmerclick.com/article/4427962654/)

# Esta función busca el mejor power:
pst <- pickSoftThreshold(datos_t, powerVector = (1:20), networkType = "unsigned") 

# Plot 1
# 
png(filename = "03_Results/01_Ajuste_Topologia.png",width = 800, height = 600, res = 120)
plot(
  pst$fitIndices[,1],
  -sign(pst$fitIndices[,3]) * pst$fitIndices[,2], #ayuda a calcular el R2
  type = "b",
  xlab = "Power",
  ylab = "Scale Free Topology Model Fit R^2",
  main = "Ajuste de Topología",
  col="blue3", pch= 16
)
abline(h = 0.8, col = "red")
# El primer punto en pasar la línea en R^2 =0.8 fue el Power = 4.
dev.off()

# Plot 2 : Conectividad media
# Validar que el power 4 retenga conexiones suficientes y no quede vacia ( conectividadpor debajo de 0)
png(filename = "03_Results/02_Conectividad_Media.png",width = 800, height = 600, res = 120)
plot(
  pst$fitIndices[, 1],
  pst$fitIndices[, 5],
  type = "b",
  xlab = "Umbral (Power)",
  ylab = "Mean Connectivity",
  main = "Conectividad Media",
  col = "darkgreen", pch = 16)
# El power 4 queda justo arriba del 0, del 5-20 estan por debajo del 0
dev.off()

# Verificación del power seleccionado
# Esta función estima el power a partir de la info obtenida con pickSoftThreshold
power <- pst$powerEstimate
power

# Se selecciona el power 4 porque es el que dadas las dos pruebas para la verificación del uso del power cumplió los requerimientos.

# Calcular la matriz de adyacencia 
# Esta permite calcular la correlación a partir de datos de expresión, se le agrega el power seleccionado para elevar los datos
# Dado que no es una red de interacción se especifica que es sin signo. 
adjacency <- adjacency(datos_t, power = 4, type = "unsigned")


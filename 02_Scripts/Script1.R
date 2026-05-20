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

# Plot 1: Ajuste de Topología:
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

# Se selecciona el power 4 porque cumple los filtros que se requieren para su selección.

# Calcular la matriz de adyacencia 
# Calcula la correlación de Pearson entre los 5,000 genes a lo largo de todas las muestras y, después, eleva cada resultado a la potencia 4. 
# Red no dirigida (no nos interesa el sentido, correlación negativa o positiva) -> sin signo -> toma el valor absoluto
# Las correlaciones débiles se desvanecen hacia el cero y las fuertes (como 0.8 o 0.9) sobreviven.
adjacency <- adjacency(datos_t, power = 4, type = "unsigned")

# Matriz de traslape topologico
# Filtro contra el ruido biologico -> Cuenta cuantos vecinos tienen en común, si comparten vecindario el enlace se vuelve más fuerte
# Compartir muchos vecinos -> solapamiento topologico alto (reduce el ruido de falsos positivos y hace que los clusters sean robustos y limpios)
TOM <- TOMsimilarity(adjacency)

# La matriz TOM perdio etiquetas entonces le ponemos los nombres de los genes de la matriz de adyacencia
rownames(TOM) <- colnames(adjacency)
colnames(TOM) <- colnames(adjacency)

# Eliminar la diagonal
# Cada gen tiene relación perfecta de 1 consigo mismo, para evitar bucles la eliminamos poniendole 0 en la diagonal
diag(TOM) <- 0

# Duplicamos la matriz TOM para filtrar
tom_filtrada <- TOM

# Filtramos la matriz para eliminar ruido de fondo, en un umbral 0.05 (arbitrario)
# Todo lo que este por debajo de este umbral se le pondra un 0
tom_filtrada[tom_filtrada < 0.05] <- 0


# Transformar datos para realizar la red compleja
red_igraph_tom <- graph_from_adjacency_matrix(tom_filtrada, # transforma la matriz plana de números en un objeto matemático de red (para trabajar con nodos y enlaces)
                                              mode = "undirected",  # los enlaces no tienen dirección (red no dirigida)
                                              weighted = TRUE, # Guarda el valor numérico de la matriz exacto de la TOM (el peso de cada linea), en lugar de tratarlo como datos simples (si/no)
                                              diag = FALSE) # ignora la diagonal de la matriz para evitar que los nodos se conecten consigo mismos

# Detección de módulos (genes que se coexpresan juntos ya que participan en la misma via biológica)
comunidades_tom <- cluster_louvain(red_igraph_tom) # agrupa los nodos altamente conectados entre si
grupos_tom <- membership(comunidades_tom) # indica a que grupo pertenece cada gen (permite clasificar)
grupos_tom


## Aislar aquellos genes con más conexiones
grados <- degree(red_igraph_tom) # cuenta cuantas conexiones tiene cada gen
grados

# seleccionar 150 de los genes que esten más conectados (hubs)
genes_hubs_tom <- names (sort(grados, decreasing = TRUE)[1:150]) # ordena los genes de mayor a menor según su número de conexiones. (names se queda únicamente con los nombres de los números)
genes_hubs_tom

# Hacer la red con esos 150 genes
sub_red_tom <- subgraph(red_igraph_tom, genes_hubs_tom) # subgraph lo que hace es extraer un subconjunto de una red más grande
# toma la red grande, busca esos 150 genes y elimina todo lo demás


# Hacer el plot de esa red con los 150 genes y alojarlos en results
png(filename = "03_Results/03_Red.png",width = 800, height = 600, res = 120)

plot (sub_red_tom,
      vertex.size = 5, # el tamaño de los nodos
      vertex.color = grupos_tom[genes_hubs_tom], # define el color de cada nodo
      vertex.label = NA, # Oculta los nombres de los genes (para que no se encimen)
      edge.width = E(sub_red_tom)$weight * 1.5, # Controla el grosor de las lineas. Extrae el valor de la fuerza de conexión de la matriz TOM (va de 0 a 1). Al multiplicarlo por 1.5 las conexiones más fuertes se ven más gruesas.
      edge.color = "grey80",  # Enlaces gris claro de fondo para estética limpia
      layout = layout_with_fr,  # Sirve para visualizar redes. Agrupa los genes relacionados en modulos
      main = "Red de Coexpresión basada en Matriz TOM (Top 150 Hubs)")
dev.off()





#Instalar paquetes a utilizar
BiocManager::install("WGCNA", force = TRUE)
library(WGCNA)
install.packages("igraph")
library(igraph)

#Investigamos y esta función ignora los demás metadatos
#datos <- read.delim("01_RawData/GSE25066_series_matrix.txt",
                    comment.char = "!",
                    header = TRUE)

# Pasamos la primera columna (Sondas/Genes) a los nombres de las filas (rownames)
#rownames(datos) <- datos[, 1]
#datos <- datos [,-1]

#Varianza 
#varianzas <- apply(datos,1,var)

#genes_variables <- names(sort(varianzas, decreasing = TRUE)[1:5000])

#datos_filtrados <- datos[genes_variables, ]


### Lo guardamos en una objeto RDS para que no tuvieran que hacer el procesamineto previo ya que ya base de datos estaba muy pesada, pero en el README esta el link 

#saveRDS(datos_filtrados, file = "01_RawData/datos.rds") #esto ya no lo tienen que hacer porque ya lo guardé en la carpeta 
datos_filtrados<- readRDS( file = "01_RawData/datos.rds" )

#transposición, para cambiar los genes de fila a columna, y las muestras de columna a fila
datos_t <- t(datos_filtrados)

# POWER: para determinar el numero por el cual se hara la potencia
# Esto fortalece la correlación fuerte y debilita la correlación débil o negativa, haciendo que el valor de correlación sea más consistente conRed sin escalaLas características son biológicamente más significativas. 
#Esta función busca el mejor power:
pst <- pickSoftThreshold(datos_t, powerVector = (1:20), networkType = "unsigned") 

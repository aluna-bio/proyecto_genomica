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


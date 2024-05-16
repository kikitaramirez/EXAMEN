#Creador: Ana Paula Ramirez Avila 
#Analisis de PCR 

#install.packages("pacman") 

library(pacman) 

p_load("readr", #para llamar a las bases de datos 
       "dplyr") #para facilitar el manejo de datos

##########################################

datos_pcr <- read_csv(file = "https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/examen2")

head(datos_pcr)

#############################################
#Obtención de los genes referencia y de interes

actina <- datos_pcr %>%
  slice(1)

actina

genes_interes <- datos_pcr %>%
  slice(-1)

genes_interes

##########################################
# promediar los controles y los tratamientos

promedio_actina <- actina %>%
  mutate(Mean_Cx=(Cx1+Cx2+Cx3)/3) %>%
  mutate(Mean_Tx=(T1+T2+T3)/3) 

promedio_actina

promedio_GI <- genes_interes %>%
  mutate(Mean_Cx=(Cx1+Cx2+Cx3)/3) %>%
  mutate(Mean_Tx=(T1+T2+T3)/3)  

promedio_GI

########################################

# Análisis DCt

DCT <- promedio_GI %>%
  mutate(DCT_Cx=(Mean_Cx-promedio_actina$Mean_Cx)) %>%
  mutate(DCT_Tx=(Mean_Tx-promedio_actina$Mean_Tx)) 

DCT
###########################

# Análisis DDCT

DDCT <- DCT %>%
  mutate(DDCT=(DCT_Tx-DCT_Cx)) %>%
  mutate("2^-DDCT"=(2^(-DDCT)))

DDCT

#######################

#Guardar tabla

write.csv(DDCT, "2DDCTFinal.csv", row.names = FALSE)


# LauraMonleon_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento
getwd()
# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)
paquetes <- c("readr", "ggplot2") # Se importan los archivos y se leen
instalar <- paquetes[!paquetes %in% installed.packages()]
if(length(instalar)) install.packages(instalar) # Pongo if por si alguno de los archivos no está instalado, que se instales y así evitar problemas y errores de código
lapply(paquetes, library, character.only = TRUE) # lapply() sube todos los paquetes para poder trabajar con ellos
datos <- read_csv("datos_biomed.csv") 

# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)
head(datos) # primeras filas para ver como son los datos

summary(datos) # resumen de las variables

dim(datos)

str(datos) # observar como son los tipos de datos

# Para saber cuántas variables (columnas) hay
ncol(datos)
# Observar cuántos tratamientos hay 
length(unique(datos$Tratamiento))
# Observar exactamente qué tratamientos son
unique(datos$Tratamiento)

# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)
# Genero boxplots para comparar las variables que se quieren observar en función de los diferentes tratamientos

par(mfrow=c(1,3)) # Lo divido en 3 para obtener 3 gráficos distintos
boxplot(Glucosa ~ Tratamiento, data = datos,
        col = c("hotpink", "purple", "lightpink"), # Elijo los colores con los que trabajar
        main = "Distribución de Glucosa por tratamiento", # Elijo el título del gráfico
        xlab = "Tratamientos", ylab = "Glucosa (mg/dL)") # Nombro los ejes X e Y con las variables

boxplot(Presion ~ Tratamiento, data = datos,
        col = c("blue", "orange", "yellow"),
        main = "Distribución de Presión por tratamiento",
        xlab = "Tratamientos", ylab = "Presión (mmHg)")

boxplot(Colesterol ~ Tratamiento, data = datos,
        col = c("red", "green", "brown"),
        main = "Distribución de Colesterol por tratamiento",
        xlab = "Tratamientos", ylab = "Colesterol (mg/dL)")
# 4. Realiza un violin plot (investiga qué es). (1 pt)
# Un violin plot es un tipo de gráfico que combina la forma de un boxplot con la forma de distribución de los datos. De esta forma se puede observar la distribución completa

install.packages("vioplot") # Instalo el vioplot y lo abro
library(vioplot)
par(mfrow = c(1, 3)) # Vuelvo a dividirlo en 3 para obtener 3 resultados

# Glucosa
vioplot(Glucosa ~ Tratamiento, data = datos,
        col = c("blue", "pink", "green"),
        main = "Distribución de Glucosa por tratamiento",
        xlab = "Tratamientos", ylab = "Glucosa (mg/dL)")

# Presión
vioplot(Presion ~ Tratamiento, data = datos,
        col = c("blue", "pink", "green"),
        main = "Distribución de Presión por tratamiento",
        xlab = "Tratamientos", ylab = "Presión (mmHg)")

# Colesterol
vioplot(Colesterol ~ Tratamiento, data = datos,
        col = c("blue", "pink", "green"),
        main = "Distribución de Colesterol por tratamiento",
        xlab = "Tratamientos", ylab = "Colesterol (mg/dL)")


# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte inferior derecha. (1 pt)
plot(datos$Glucosa, datos$Presion, # Plot() se usa para generar gráficos sencillos
     col = as.numeric(factor(datos$Tratamiento)),  # Se define un color para cada tratamiento
     pch = 19,    # con pch se generan puntos en el gráfico para estudiar así la dispersión. pongo 19 porque es el tamaño más adecuado                                   
     main = "Glucosa vs Presión",
     xlab = "Glucosa (mg/dL)",
     ylab = "Presión (mmHg)")

legend("bottomright", # Función legend() para poner los nombres de los tratamientos
       legend = levels(factor(datos$Tratamiento)),  
       col = 1:length(unique(datos$Tratamiento)),  
       pch = 19, 
       title = "Tratamientos")

# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)
# Facet Grid sirve para hacer varios gráficos a la vez separándolos por alguna variable y los organiza en cuadrículas
library(ggplot2) # Repito el gráfico anterior pero con ggplot para poder personalizarlo más
ggplot(datos, aes(x = Presion, y = Colesterol, color = Tratamiento)) +
  geom_point(size = 2) +   # geom_point() Dibuja puntos que representan valores y su relación entre ellas                 
  geom_smooth(method = "lm", se = FALSE) +   # geom_smooth() Para añadir una linea de tendencia
  facet_wrap(~ Tratamiento) +    # facet wrap sirve para definir la variable y las dispone automáticamente            
  labs(title = "Relación entre Colesterol y Presión por tratamiento", #labs() se pone al principio para después personalizar el gráfico (el título, colores, los ejes, etc.)
       x = "Presión (mmHg)",
       y = "Colesterol (mg/dL)") +
  theme_minimal() # theme_minimal() para hacerlo más simplificado

# 7. Realiza un histogramas para cada variable. (0.5 pts)
par(mfrow = c(1, 3)) # De nuevo, se divide en tres
hist(datos$Glucosa,
     main = "Histograma de Glucosa",
     col = "blue",
     xlab = "Glucosa (mg/dL)")
hist(datos$Presion,
     main = "Histograma de Presión",
     col = "pink",
     xlab = "Presión (mmHg)")
hist(datos$Colesterol,
     main = "Histograma de Colesterol",
     col = "green",
     xlab = "Colesterol (mg/dL)")

# 8. Crea un factor a partir del tratamiento. Investifa factor(). (1 pt)
factor(datos$Tratamiento) # Convierto la columna de Tratamiento en factor para poder usarlo en los análisis y gráficos sin errores
levels(factor(datos$Tratamiento)) # Se usa levels() para ver los nombres de los niveles que tiene el factor

# 9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5 pts)
aggregate(Glucosa ~ Tratamiento, data = datos, FUN = mean) # Con mean calculo la media
aggregate(Glucosa ~ Tratamiento, data = datos, FUN = sd) # con sd calculo la desviación estándar
# aggregate() sirve para ejecutar una función (ej: mean/sd)

# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)
# Guardo los datos de cada tratamiento por separado para poder analizarlos mejor y de forma más cómoda
placebo  <- subset(datos, Tratamiento == "Placebo")
farmacoA <- subset(datos, Tratamiento == "FarmacoA")
farmacoB <- subset(datos, Tratamiento == "FarmacoB")
# Uso subset() para extraer solo las filas que correspondan a cada uno de los tratamientos

head(placebo) # con head() se pueden ver las primeras líneas del código, en este caso del tratamiento con placebo

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)
# Los datos son normales (siguen una distribución normal) ya que su p-value es p >0,05
# Como los datos son normales, uso t-test para comparar medias entre tratamientos
# Lo hago para glucosa y para presión y colesterol, ya que de esta manera se puede comprobar si los tratamientos también afectan a estas variables
t.test(placebo$Glucosa,farmacoA$Glucosa)
t.test(placebo$Glucosa,farmacoB$Glucosa)
t.test(farmacoB$Glucosa,farmacoA$Glucosa)

# Uso el test de Shapiro-Wilk para ver si la variable Glucosa sigue una distribución normal, y como su p-value es p>0,05, asumo que los datos son normales
shapiro.test(placebo$Presión) # En esta línea de código y las siguietes, se selecciona únicamente una variable o columna en concreto
shapiro.test(farmacoA$Presión)
shapiro.test(farmacoB$Presión)

t.test(placebo$Glucosa,farmacoA$Presión)
t.test(placebo$Glucosa,farmacoB$Presión)
t.test(farmacoB$Glucosa,farmacoA$Presión)

shapiro.test(placebo$Colesterol)
shapiro.test(farmacoA$Colesterol)
shapiro.test(farmacoB$Colesterol)

t.test(placebo$Glucosa,farmacoA$Colesterol)
t.test(placebo$Glucosa,farmacoB$Colesterol)
t.test(farmacoB$Glucosa,farmacoA$Colesterol)

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)
# uso ANOVA para comparar las medias de los tres tratamientos al mismo tiempo
anova_glucosa <- aov(Glucosa ~ Tratamiento, data = datos) # aov sirve para detectar diferencias globales entre los grupos
summary(anova_glucosa) # con summary se obtiene un resumen rápido con lo más importante y de forma ordenada
TukeyHSD(anova_glucosa) # 	Uso TurkeyHSD() para saber entre qué tratamientos están las diferencias de los grupos



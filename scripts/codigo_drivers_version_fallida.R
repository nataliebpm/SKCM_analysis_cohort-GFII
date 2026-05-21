# Cargar el dataset de mutaciones
# Este Dataset ya se encuentra anotado, filtrado y el variant calling ya fue realizado.

# Cargar librerias
library(maftools)
library(GenomeInfoDb)
library(dndscv)

# Cargar el MAF y definir la ruta donde se encuentra el archivo
SKCM <- read.maf(maf = "/Users/monicareyes/Documents/Semestre4_LCG/Genomica_funcional2/Tareas_proyectos/Proyecto_final/SKCM/data_mutations.txt")
muts <- SKCM@data

# Filtramos los datos que nos van a servir para hacer el análisis con dndscv
muts_dndscv <- data.frame(
  sampleID = muts$Tumor_Sample_Barcode,
  chr      = muts$Chromosome,
  pos      = muts$Start_Position,
  ref      = muts$Reference_Allele,
  alt      = muts$Tumor_Seq_Allele2
)

# Visualizamos como quedaron los datos
head(muts_dndscv)

# Numero de muestras unicas que tenemos
length(unique(muts_dndscv$sampleID))

# Numero de mutaciones registradas
nrow(muts_dndscv)

# Revisamos si hay hipermutadores (>500 mutaciones)
barplot(sort(table(muts_dndscv$sampleID)), 
        ylab = "Number of mutations",
        xlab = "Donors", 
        las = 2, 
        names.arg = "")
abline(h = 500, col = "red", lty = 2)

# Contar mutaciones por muestra
muts_por_muestra <- table(muts_dndscv$sampleID)

# Identificar hipermutadores (>500 mutaciones)
hipermutadores <- muts_por_muestra[muts_por_muestra > 500]

# Cuántos hay
length(hipermutadores)

# dndscv excluye automáticamente los hipermutadores
dout = dndscv(muts_dndscv,
              max_muts_per_gene_per_sample = 3,
              max_coding_muts_per_sample = 500,
              outmats = T)
names(dout)

# Tabla de genes significativos 
dout$sel_cv[which(dout$sel_cv$qglobal_cv < 0.1),]

# Visualizar los oncogenes y supresores tumorales encontrados
drivers <- dout$sel_cv[which(dout$sel_cv$qglobal_cv < 0.1),]
drivers$type <- ifelse((drivers$n_non + drivers$n_spl + drivers$n_ind) / 
                         (drivers$n_mis + drivers$n_non + drivers$n_spl + drivers$n_ind) > 0.20, 
                       "Tumor Suppressor", "Oncogene")
oncogenes <- drivers$gene_name[drivers$type == "Oncogene"]
supresores <- drivers$gene_name[drivers$type == "Tumor Suppressor"]

max_len <- max(length(oncogenes), length(supresores))
tabla <- data.frame(
  Oncogenes = c(oncogenes, rep(NA, max_len - length(oncogenes))),
  Tumor_Suppressors = c(supresores, rep(NA, max_len - length(supresores)))
)
print(tabla)

# Selección positiva y negativa basada en wmis_cv
drivers$seleccion <- ifelse(drivers$wmis_cv > 1, "Selección positiva", "Selección negativa")

tabla_seleccion <- data.frame(
  Gen        = drivers$gene_name,
  Tipo       = drivers$type,
  Seleccion  = drivers$seleccion,
  wmis_cv    = round(drivers$wmis_cv, 3),
  wnon_cv    = round(drivers$wnon_cv, 3),
  qglobal_cv = format(drivers$qglobal_cv, scientific = TRUE, digits = 3)
)
tabla_seleccion <- tabla_seleccion[order(-drivers$wmis_cv), ]
print(tabla_seleccion)

# Señales globales de selección
print(dout$globaldnds)

# Análisis de hotspots
dout$annotmuts$gene_and_aachange = paste(dout$annotmuts$gene,
                                         dout$annotmuts$aachange,
                                         dout$annotmuts$ntchange,
                                         dout$annotmuts$pos,
                                         dout$annotmuts$impact, sep=":")
sort(table(dout$annotmuts$gene_and_aachange), decreasing = T)[1:10]

# sitednds - solo mutaciones en Cancer Gene Census genes (v81)
data("cancergenes_cgc81", package = "dndscv")
dout_cancergenes = dndscv(muts_dndscv, outmats = T, gene_list = known_cancergenes)
sout = sitednds(dout_cancergenes)
sout$recursites[which(sout$recursites$qval < 0.1),]

# Análisis de selección a nivel de codón
data("refcds_hg19", package = "dndscv")
RefCDS_codon = buildcodon(RefCDS)
codon_dnds = codondnds(dout_cancergenes, 
                       RefCDS_codon,
                       theta_option = "conservative", 
                       min_recurr = 2)
codon_dnds$recurcodons[which(codon_dnds$recurcodons$qval < 0.1),]

muts_drivers <- dout$annotmuts[dout$annotmuts$gene %in% drivers$gene_name, ]

# Contar pacientes únicos por gen
pacientes_por_gen <- tapply(muts_drivers$sampleID, 
                            muts_drivers$gene, 
                            function(x) length(unique(x)))

# Crear la tabla
tabla_drivers <- data.frame(
  Gen              = names(pacientes_por_gen),
  Numero_pacientes = as.numeric(pacientes_por_gen)
)

# Ordenar de mayor a menor frecuencia
tabla_drivers <- tabla_drivers[order(-tabla_drivers$Numero_pacientes), ]

# Verificar que se creó bien
cat("tabla_drivers creada con", nrow(tabla_drivers), "genes\n")
print(tabla_drivers)

##### Filtrado de genes ########

# En esta sección vamos a filtrar unicamente aquellos genes que tengan una selección positiva estadistica
#(ya que niguno de estos genes esta bajo selección negativa), que tenga hotspots recurrentes y que tengan una
# frecuencia alta en la cohorte, de este modo los genes que cumplan los tres filtros, serán significativos para 
# nuestro trabajo


# Filtramos primeramente los genes que son muy largos, por que adquieren muchas mutaciones por ser muy grandes en pares de bases
# y no por se significativos realmente

genes_excluir <- c("TTN", "MUC16", "DNAH5")


# evidencia 1: Selección positiva estadística (dout$sel_cv)

evidencia1 <- dout$sel_cv[dout$sel_cv$qglobal_cv < 0.1, ]
evidencia1 <- evidencia1[!evidencia1$gene_name %in% genes_excluir, ]
genes_ev1  <- evidencia1$gene_name

cat("Genes con selección positiva (qglobal_cv < 0.1):\n")
print(genes_ev1)


# evidencia 2: Hotspots recurrentes (sout$recursites)

evidencia2 <- sout$recursites[sout$recursites$qval < 0.1, ]
evidencia2 <- evidencia2[!evidencia2$gene %in% genes_excluir, ]
genes_ev2  <- unique(evidencia2$gene)

cat("\nGenes con hotspots recurrentes (qval < 0.1):\n")
print(genes_ev2)


# evidencia 3: Frecuencia en la cohorte (tabla_drivers_filtrada)
#
# Usamos la mediana de pacientes como umbral mínimo de frecuencia, por que la media es muy sensible a outliers.
umbral_freq <- median(tabla_drivers_filtrada$Numero_pacientes)
evidencia3  <- tabla_drivers_filtrada[tabla_drivers_filtrada$Numero_pacientes >= umbral_freq, ]
evidencia3  <- evidencia3[!evidencia3$Gen %in% genes_excluir, ]
genes_ev3   <- evidencia3$Gen

cat("\nGenes frecuentes en la cohorte (>= mediana de pacientes):\n")
print(genes_ev3)
cat("Umbral usado:", umbral_freq, "pacientes\n")


#### Resultados ###

# Genes en que cumplen con las 3 evidencias — drivers más robustos
drivers_robustos <- Reduce(intersect, list(genes_ev1, genes_ev2, genes_ev3))

#Genes en al menos 2 evidencias
drivers_2de3 <- unique(c(
  intersect(genes_ev1, genes_ev2),
  intersect(genes_ev1, genes_ev3),
  intersect(genes_ev2, genes_ev3)
))
drivers_2de3 <- drivers_2de3[!drivers_2de3 %in% drivers_robustos]  # quitar los que ya están en 3/3

print(drivers_2de3)
cat("drivers con las 3 evidencias (más robustos):\n")
print(drivers_robustos)


# Resumen final
todos_los_genes <- unique(c(genes_ev1, genes_ev2, genes_ev3))

tabla_final <- data.frame(
  Gen            = todos_los_genes,
  Sel_positiva   = todos_los_genes %in% genes_ev1,   # Evidencia 1
  Hotspot        = todos_los_genes %in% genes_ev2,   # Evidencia 2
  Frecuente      = todos_los_genes %in% genes_ev3,   # Evidencia 3
  N_evidencias   = rowSums(data.frame(
    todos_los_genes %in% genes_ev1,
    todos_los_genes %in% genes_ev2,
    todos_los_genes %in% genes_ev3
  ))
)

# Agregar tipo (oncogén o supresor) desde la tabla drivers
tabla_final$Tipo <- drivers$type[match(tabla_final$Gen, drivers$gene_name)]

# Agregar número de pacientes
tabla_final$N_pacientes <- tabla_drivers_filtrada$Numero_pacientes[
  match(tabla_final$Gen, tabla_drivers_filtrada$Gen) 
]

# Ordenar por número de evidencias y luego por pacientes
tabla_final <- tabla_final[order(-tabla_final$N_evidencias, -tabla_final$N_pacientes), ]

print(tabla_final)

###### Oncoplot solamente con los drivers de alta confianza #######

oncoplot(maf = SKCM,
         genes = drivers_robustos,
         genesToIgnore = c("TTN", "MUC16", "DNAH5"),
         draw_titv = TRUE)

##### Lolliplot de los drivers ######
drivers_robustos <- c("BRAF", "NRAS", "B2M", "KIT", "GNAQ", "MAP2K1", "GNA11", "CDK4")

for (gen in genes_drivers_robustos) {
  
  # Verificar que el gen tiene suficientes mutaciones para graficar
  n_muts <- sum(SKCM@data$Hugo_Symbol == gen)
  
  if (n_muts < 2) {
    cat("Gen", gen, "tiene menos de 2 mutaciones, saltando...\n")
    next
  }
  
  cat("Graficando:", gen, "(", n_muts, "mutaciones )\n")
  
  lollipopPlot(
    maf            = SKCM,
    gene           = gen,
    showMutationRate = TRUE,   
    labelPos       = "all",    # se eiqueta todas las posiciones
    labPosAngle    = 45        # ángulo de las etiquetas
  )
}
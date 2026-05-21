# Cargar el dataset de mutaciones
# Este Dataset ya se encuentra anotado, filtrado y el variant calling ya fue realizado.

# Cargar librerias
library(maftools)
library(GenomeInfoDb)
library(dndscv)

# Cargar el MAF y definir la ruta donde se encuentra el archivo
SKCM <- read.maf(maf = "/Users/monicareyes/Documents/Semestre4_LCG/Genomica_funcional2/Tareas_proyectos/Proyecto_final/SKCM/data_mutations.txt")

# Combinar mutaciones codificantes Y silenciosas
# (dndscv necesita ambas para calcular correctamente dN/dS)
todas <- rbind(SKCM@data, SKCM@maf.silent)  # <-- línea que faltaba

# Filtramos los datos que nos van a servir para hacer el análisis con dndscv
muts_dndscv <- data.frame(
  sampleID = todas$Tumor_Sample_Barcode,
  chr      = todas$Chromosome,
  pos      = todas$Start_Position,
  ref      = todas$Reference_Allele,
  alt      = todas$Tumor_Seq_Allele2
)

# Filtrar solo SNVs
muts_dndscv <- muts_dndscv[
  nchar(muts_dndscv$ref) == 1 & nchar(muts_dndscv$alt) == 1, 
]

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
dout$sel_cv[which(dout$sel_cv$qglobal_cv < 0.2),]

# Tabla de genes significativos
genes_sig <- dout$sel_cv[which(dout$sel_cv$qglobal_cv < 0.2), ]

# Extraer solo los nombres
drivers <- genes_sig$gene_name

###### Oncoplot solamente con los drivers de alta confianza #######
oncoplot(maf          = SKCM,
         genes        = drivers,
         genesToIgnore = c("TTN", "MUC16", "DNAH5"),
         draw_titv    = TRUE)

##### Lolliplot de los drivers ######
for (gen in drivers) {
  
  n_muts <- sum(SKCM@data$Hugo_Symbol == gen)
  
  if (n_muts < 2) {
    cat("Gen", gen, "tiene menos de 2 mutaciones, saltando...\n")
    next
  }
  
  cat("Graficando:", gen, "(", n_muts, "mutaciones)\n")
  
  lollipopPlot(
    maf              = SKCM,
    gene             = gen,
    showMutationRate = TRUE,
    labelPos         = "all",
    labPosAngle      = 45
  )
}
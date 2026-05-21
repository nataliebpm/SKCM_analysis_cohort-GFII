###########################################################################################################################################################
# date: 12 mayo 2026
# Proyecto de Genómica Funcional II
# Análisis de datos mutacionales para melanoma

rm(list = ls())

#------------------------------------------1. Reformatear el archivo de mutaciones---------------------------------------------------------------------

library(tidyverse)

# Leer tu archivo de mutaciones de melanoma
maf_cbioportal = read.delim('/home/user/Ciencias_Genómicas_LCGEJ/Genómica_Funcional/Genomica_Funcional_II/Proyecto_Final/data_mutations.txt')

# Seleccionar solo las columnas que necesita SigProfiler
maf_sp = maf_cbioportal %>%
  select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome,
         Start_Position, End_Position, Strand, Variant_Classification,
         Variant_Type, Reference_Allele, Tumor_Seq_Allele1,
         Tumor_Seq_Allele2, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode)

# Filtrar solo sustituciones de una base (SNPs) — descartar indels
maf_sp = maf_sp %>%
  filter(Variant_Type == 'SNP')

# Crear carpetas para los resultados
dir.create('signatures')
dir.create('signatures/SPMG/')

# Guardar el archivo reformateado
write.table(maf_sp, 'signatures/SPMG/data_mutations.maf',
            quote = F, row.names = F, sep = '\t')

#----------------------------------2. Generar la matriz mutacional SBS96------------------------------------------------------------------------

# 1. Definimos la ruta ANTES de cargar nada
#Sys.setenv(RETICULATE_PYTHON = "/home/user/miniconda3/envs/mutational_signatures/bin/python")

# 2. Ahora cargamos las librerías
library(reticulate)
#py_config()
library(SigProfilerMatrixGeneratorR)

# 3. Forzamos la inicialización
use_condaenv("/home/user/miniconda3/envs/mutational_signatures", required = TRUE)

# 4. Verificación (Debe salir la ruta de tu miniconda)
#py_config()

# 4. Verificación rápida (opcional pero recomendada)
# Si esto falla, el problema es la instalación de Python, no el código de R
#if(!py_module_available("SigProfilerMatrixGenerator")){
  #stop("El módulo de Python no está en el entorno. Revisa: pip install SigProfilerMatrixGenerator")
#}

# 5. Generar matrices
# Nota: La primera vez descargará el genoma GRCh37 si no está presente
matrices <- SigProfilerMatrixGeneratorR(
  project = "SKCM", 
  genome  = "GRCh37",
  matrix_path = "./signatures/SPMG",
  plot = FALSE,
  exome = TRUE
)

#----------------------------------3. Visualizar perfiles mutacionales------------------------------------------------------------------------

library(SigProfilerPlottingR)

# Perfil de todas las muestras
plotSBS(
  matrix_path = 'signatures/SPMG/output/SBS/SKCM.SBS96.exome',
  output_path = 'signatures/SPMG/output/SBS/',
  project     = 'SKCM',
  plot_type   = '96',
  percentage  = FALSE
)

#----------------------------------4. Perfil promedio en la cohorte------------------------------------------------------------------------

mut_matrix = matrices[['96']]

# Convertir a valores relativos (para no sesgar por muestras muy mutadas)
relative_mut_matrix = apply(mut_matrix, 2, prop.table)

# Calcular promedio
average_mut_matrix = rowMeans(relative_mut_matrix)
average_mut_matrix = data.frame(Average_SKCM = average_mut_matrix)

# Guardar y graficar
average_mut_matrix_to_print = cbind(rownames(average_mut_matrix), average_mut_matrix)
colnames(average_mut_matrix_to_print)[1] = 'MutationType'

write.table(average_mut_matrix_to_print, 'signatures/avg_SKCM.SBS96.all',
            quote = F, row.names = F, sep = '\t')

plotSBS(
  matrix_path = 'signatures/avg_SKCM.SBS96.all',
  output_path = 'signatures/',
  project     = 'avg_SKCM',
  plot_type   = '96',
  percentage  = TRUE
)

#----------------------------------4. Perfil promedio en subgrupos------------------------------------------------------------------------

# Leer el archivo clínico
#metadata = read.delim('/home/user/Ciencias_Genómicas_LCGEJ/Genómica_Funcional/Genomica_Funcional_II/Proyecto_Final/skcm_tcga_pan_can_atlas_2018_clinical_data.tsv')

# Filtrar para usar solo pacientes con datos de mutación
#metadata = metadata %>%
 # filter(Sample.ID %in% maf_sp$Tumor_Sample_Barcode)

# Ver cuántas muestras hay por subgrupo
#table(metadata$Subtype)


#----------------------------------4. Extraer firmas de novo------------------------------------------------------------------------

library(SigProfilerExtractorR)

# Aseguramos que reticulate use el entorno correcto
reticulate::use_condaenv("/home/user/miniconda3/envs/mutational_signatures", required = TRUE)
sigprofilerextractor(
  input_type         = 'matrix',
  output             = 'signatures/SPE/',
  input_data         = 'signatures/SPMG/output/SBS/SKCM.SBS96.exome',
  nmf_replicates     = 3,    # en producción real usar 100
  minimum_signatures = 1,
  maximum_signatures = 10,   # en melanoma pueden haber varias firmas activas
  exome              = T
)


#----------------------------------4. Asignar firmas COSMIC a cada muestra------------------------------------------------------------------------

library(SigProfilerAssignmentR)

cosmic_fit(
  samples    = 'signatures/SPMG/output/SBS/SKCM.SBS96.exome',
  output     = 'signatures/SPA',
  input_type = 'matrix',
  exome      = T
)


#----------------------------------4. Resultados------------------------------------------------------------------------

# Leer actividades de firmas (resultado de SigProfilerExtractor)
acts = read.delim('signatures/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

# Promedio de mutaciones por firma en todo el cohorte
avg_acts = colMeans(acts[, -1])
barplot(avg_acts, las=2, main="Actividad promedio de firmas - SKCM")

# Reformatear para ggplot
acts_tidy = acts %>%
  pivot_longer(cols = !Samples, names_to = 'Signature', values_to = 'Mutations')

# Unir con metadatos
acts_and_metadata = acts_tidy %>%
  rename(Sample.ID = Samples) %>%
  left_join(metadata)

# Promedio de actividad por subgrupo
acts_per_subgroup = acts_and_metadata %>%
  group_by(Subtype, Signature) %>%
  summarise(Avg_mutations = mean(Mutations))

# Graficar por subgrupo
ggplot(acts_per_subgroup) +
  aes(x = reorder(Subtype, Avg_mutations), y = Avg_mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  labs(x = 'Subtipo de Melanoma', y = 'Proporción de mutaciones')

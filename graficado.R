library("dplyr")
library("ggplot2")
library("tidyr")

# leer la data
counts <- read.csv( file = "all.counts.tsv", sep = "\t",
                    header = T, skip = 1 )

# Obtener nombre de genes ordenados
ordenados <- counts$Geneid

# Quedarse solo con columnas utiles
selected <- counts %>% 
  select( -Chr, -Start, -End, -Strand )

# lo pasamos a formato largo
long <- gather( selected,
                key = "Sample",
                value = "Counts",
                -Geneid, -Length )

## leemos los datos sobre totales de lecturas
totales <- read.csv( file = "all.counts.tsv.summary", sep = "\t",
                     header = T )

# lo pasamos a formato largo
totales_long <- gather( totales,
                        key = "Sample",
                        value = "Reads",
                        -Status ) %>% 
  group_by( Sample ) %>% 
  summarise( Totals = sum(Reads) ) %>% 
  ungroup( )

# agregamos los datos de lecturas totales a los conteos
unidos <- left_join( x = long,
                     y = totales_long, 
                     by = "Sample" )

# normalizamos por RPKM y por TPM
# Usamos: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
finales <- unidos %>% 
  mutate( Total_millions = Totals / 1e6,
          RPM = Counts / Total_millions,
          Length_kb = Length/ 1000,
          RPKM = RPM / Length_kb ) %>% 
  group_by( Sample ) %>% 
  mutate( RPK = Counts / Length_kb,
          Total_RPK_million = sum(RPK) / 1e6,
          TPM = RPK / Total_RPK_million )

# Creamos una funcion para graficar las lineaspor gen

plot_lines <- function( la_data, la_columna ) {

ggplot( data = la_data,
        mapping = aes( x = Geneid,
                       y = ycolumn,
                       group = Sample,
                       color = Sample ) ) +
  geom_line( size = 1 ) +
  geom_point( shape = 21, size = 2, fill = "white" ) +
  scale_x_discrete( limits = ordenados ) +
  labs( title = "Detección de Genes SARS-CoV-2 por RNA-Seq",
        subtitle = "Celulas Calu-3, epitelio de pulmón",
        x = "Gen viral",
        y = la_columna ) +
  theme_classic( base_size = 30 ) +
  theme( axis.text.x = element_text( angle = 90 ) )
  
  # save the plot
  ggsave( filename = paste0(la_columna, ".png"),
          width = 14,
          height = 10 )
  
}

# Plot counts
rename( finales,
        ycolumn = Counts ) %>% 
plot_lines( la_data = ., la_columna = "Counts" )

# Plot RPKM
rename( finales,
        ycolumn = RPKM ) %>% 
  plot_lines( la_data = ., la_columna = "RPKM" )

# Plot TPM
rename( finales,
        ycolumn = TPM ) %>% 
  plot_lines( la_data = ., la_columna = "TPM" )

# Guardar la tabla de conteos
write.csv( x = finales, file = "all.counts_with_RPKM_and_TPM.tsv", quote = T, row.names = F)

# Fin del Script :)
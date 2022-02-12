#!/bin/bash

# Hagamos el QC control de calidad de la data original
fastqc -o . inputdata/*.fastq.gz

# Hagamos el Corte o Trimming de los datos
# Para eliminar secuencias feas
trimmomatic SE \
        -phred33 \
        inputdata/infected.fastq.gz \
        infected.trim.fastq.gz \
        HEADCROP:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36 \
	2> infected.trim.log

# borramos la prueba
rm *.trim.*

# Bien, eso fue para una sola muestra
# Construyamos un loop for, que haga el Trimming para todos los fasq.gz en inputdata/*.fastq.gz
for ifile in inputdata/*.fastq.gz
do
	echo "..Trimming $ifile"
	trimmomatic SE \
        -phred33 \
        $ifile \
        $ifile.trim.fastq.gz \
        HEADCROP:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36 \
	2> $ifile.trim.log
	mv $ifile.trim.*   .
done

# Volvemos a correr QC pero sobre los datos trimmeados
fastqc -o . *.trim.fastq.gz

# Descargamos el genoma de referencia viral y construimos el indice (el mapa base)

wget http://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz

gunzip -f Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz

subread-buildindex \
        -o sarscov \
        Sars_cov_2.ASM985889v3.dna.toplevel.fa

# Ahora Mapeamos (alineamos) una muestra
subread-align \
        -t 0 \
        -i sarscov \
        -r infected.fastq.gz.trim.fastq.gz \
        --SAMoutput \
        -o infected.map.sam

# Borramos la prueba
rm *.sam*

# Bien, eso fue para una sola muestra
# Construyamos un loop for, que haga el Mapeo para todos los .trim.fastq.gz en *.trim.fastq.gz

for ifile in *.trim.fastq.gz
do
	echo "..Mapping $ifile"
	subread-align \
        -t 0 \
        -i sarscov \
        -r $ifile \
        --SAMoutput \
        -o $ifile.map.sam
done

# Descargamos la referencia con las coordenadas de cada gen
wget http://ftp.ensemblgenomes.org/pub/viruses/gtf/sars_cov_2/Sars_cov_2.ASM985889v3.101.gtf.gz

gunzip -f Sars_cov_2.ASM985889v3.101.gtf.gz


# utilizamos subreads para contar cuantas lecturas caen en cada gen
featureCounts \
        -g gene_name \
        -t CDS \
        -O \
        -s 1 \
        -a Sars_cov_2.ASM985889v3.101.gtf \
        -o all.counts.tsv \
	*.sam

# Graficamos
Rscript graficado.R

# creamos una carpeta para los resultados
mkdir -p results/

# Movemos los resultados
mv *.tsv* *.png results/

# borramos todos los temporales locales
rm *fastq*

# Sacamos control de calidad de los pasos
multiqc .

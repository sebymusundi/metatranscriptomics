#! /usr/bin/nextflow 

params.reads='/home/bioinfo/wastewater_analysis/raw_data/*_{R1,R2}_001.fastq'
params.outdir="/home/bioinfo/wastewater_analysis/results/"
params.reference='/home/bioinfo/wastewater_analysis/raw_data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz'

println "reads: $params.reads"
println "reference: $params.reference"
println "outdir: $params.outdir"



log.info """\
                        whole genome sequencing SARS-CoV-2
                        ------------------------------------------
                        reads:       "$params.reads"
                        reference:   "$params.reference"
                        outdir:      "$params.outdir"
                        ------------------------------------------
                        """
                .stripIndent()


// Quality control using fastqc

process fastqc {
                container 'biocontainers/fastqc:v0.11.9_cv8'
                publishDir "${params.outdir}", mode: 'copy'

                input:
                tuple val(sample_id), path(my_reads_1)

                output:
                tuple val(sample_id), 
                path("quality_control/${my_reads_1[0].baseName}_fastqc.html"), 
                path("quality_control/${my_reads_1[1].baseName}_fastqc.html"), 
                path("quality_control/${my_reads_1[0].baseName}_fastqc.zip"), 
                path("quality_control/${my_reads_1[1].baseName}_fastqc.zip")


                script:
                """
                mkdir  quality_control
                fastqc ${my_reads_1[0]} ${my_reads_1[1]}    -o quality_control

                """

}

// Removal of adapters using fastp

process fastp {
                    container 'nanozoo/fastp:latest'
                    publishDir "${params.outdir}", mode: 'copy'

                    input:
                    tuple val(sample_id), path(my_reads_1)

                    output:
                    tuple val(sample_id),  path("fastp_output/${my_reads_1[0].baseName}.trimmed.fastq"), path("fastp_output/${my_reads_1[1].baseName}.trimmed.fastq")
                    
                  

                    script:
                    """
                    mkdir -p fastp_output 
                    
                    fastp -i ${my_reads_1[0]} -I ${my_reads_1[1]} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -Q 20 \
                    -o fastp_output/${my_reads_1[0].baseName}.trimmed.fastq -O fastp_output/${my_reads_1[1].baseName}.trimmed.fastq 
                
                    """
}

//  Trim reads with a quality score less than 20 using trimmomatic

process trimmomatic {
                    container 'staphb/trimmomatic:latest'
                    publishDir "${params.outdir}", mode: 'copy'

                    input:
                    tuple val(sample_id), path(my_reads_1), path(my_reads_2)

                    output:
                   
                    tuple val(sample_id), path("trimmomatic_output/${my_reads_1.baseName}.paired.fastq"), 
                    path("trimmomatic_output/${my_reads_1.baseName}.unpaired.fastq"),
                    path("trimmomatic_output/${my_reads_2.baseName}.paired.fastq"),
                    path("trimmomatic_output/${my_reads_2.baseName}.unpaired.fastq")
                  

                    script:
                    """
                    mkdir -p trimmomatic_output 
                    
                    trimmomatic-0.39.jar PE ${my_reads_1} ${my_reads_2}  trimmomatic_output/${my_reads_1.baseName}.paired.fastq trimmomatic_output/${my_reads_1.baseName}.unpaired.fastq \
                    trimmomatic_output/${my_reads_2.baseName}.paired.fastq trimmomatic_output/${my_reads_2.baseName}.unpaired.fastq SLIDINGWINDOW:4:20 
                    """
}

// Checking the quality of trimmed files 

process fastqc_trimmed{
                        container 'biocontainers/fastqc:v0.11.9_cv8'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        tuple val(sample_id), path(trimmed_reads_1), path(trimmed_reads_2), path(trimmed_reads_3), path(trimmed_reads_4)

                        output:
                        tuple val(sample_id), 
                        path("quality_control_trimmed/${trimmed_reads_1.baseName}_fastqc.html"), 
                        path("quality_control_trimmed/${trimmed_reads_3.baseName}_fastqc.html"), 
                        path("quality_control_trimmed/${trimmed_reads_1.baseName}_fastqc.zip"), 
                        path("quality_control_trimmed/${trimmed_reads_3.baseName}_fastqc.zip")


                        script:
                        """
                        mkdir  quality_control_trimmed
                        fastqc ${trimmed_reads_1}  ${trimmed_reads_3}    -o quality_control_trimmed

                        """
}

// removal of human contaminating reads revealed by kraken  

// Creating an index for the reference genome 

process bowtie_human {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        path reference

                        output:
                        tuple val("bowtie_human"), path("bowtie_human*")
                
                        script:
                        """
                        bowtie2-build ${reference}  bowtie_human
                        """
}

// map the indexed genome against the reads 

process bowtie_align_human {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'
                        memory '16 GB'

                        input:
                        tuple val(x), path(paired_reads_1), path(unpaired_reads_1), path(paired_reads_2), path(unpaired_reads_2), val(index), path(align_input)
                    
                    
                        output:
                        tuple val(x), path("bowtie_output_human/${x}.sam")
                        
                        script:
                        """
                        mkdir -p bowtie_output_human
                        
                 
                        bowtie2 -x ${index} -1 ${paired_reads_1} -2 ${paired_reads_2} -S bowtie_output_human/${x}.sam                   
                        """
}


// Converting sam file to become Bam file followed by sorting and indexing the associated BAM file

process samtools {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files/${sam_file.baseName}.sorted.bam"), path("bam_files/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                     mkdir bam_files 
                    samtools view -O BAM ${sam_file} -o  bam_files/${sam_file.baseName}.bam
                    samtools sort  bam_files/${sam_file.baseName}.bam -o  bam_files/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files/${sam_file.baseName}.sorted.bam
                     """
}

// filter out unmapped reads 
process filter_unmapped_reads {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(bam_file), path(bam_file_index)

                     output:
                     tuple val(x), path("bam_files_unmapped/${bam_file.baseName}.unmapped.bam")

                     script:
                     """
                     mkdir bam_files_unmapped
                     samtools view -b -f 12 -F 256 ${bam_file} -o bam_files_unmapped/${bam_file.baseName}.unmapped.bam
                    
                     """
}

// split paired end reds into separated fastq files 

process split_file_sort  {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(bam_file)

                     output:
                     tuple val(x), path("bam_files_sort/${bam_file.baseName}.unmapped.sorted.bam")

                     script:
                     """
                     mkdir bam_files_sort
                     samtools sort ${bam_file} -o bam_files_sort/${bam_file.baseName}.unmapped.sorted.bam -O BAM
                           
                     """
}

// split files to fastq 

process split_file  {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(bam_file)

                     output:
                     tuple val(x), path("bam_files_split/${x}_R1_001.fastq"), path("bam_files_split/${x}_R2_001.fastq")

                     script:
                     """
                     mkdir bam_files_split
                     samtools fastq -1 bam_files_split/${x}_R1_001.fastq -2 bam_files_split/${x}_R2_001.fastq ${bam_file}
                           
                     """
}

// combine forward  and reverse reads 

process combine_reads  {
                     
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x), path(trimmed_reads_1), path(trimmed_reads_2)


                     output:
                     tuple val(x), path("combined_reads/merged_1.fastq"), path("combined_reads/merged_2.fastq")

                     script:
                     """
                     mkdir -p combined_reads
               
                     cat ${trimmed_reads_1}    > combined_reads/merged_1.fastq
                     cat ${trimmed_reads_2}    > combined_reads/merged_2.fastq
                     """
}

// Running de novo assembly for the samples using   SPADES

process spades {
                                        container 'staphb/spades:latest'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        memory '16 GB'
                                        

                                        input:
                                        tuple val(x), path(trimmed_reads_1), path(trimmed_reads_2)

                                        output:
                                        tuple val(x), path("*")


                                        script:
                                        """
                                        mkdir spades_output
                                        spades.py --meta -1 ${trimmed_reads_1} -2 ${trimmed_reads_2} -o spades_output/${x} --tmp-dir . --threads 1
                                        """
}


process bioawk {
                                        container 'lbmc/bioawk:1.0'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        
                                        

                                        input:
                                        tuple val(x), path(contig_file)

                                        output:
                                        tuple val(x), path("filtered_contigs/${x}.fasta")


                                        script:
                                        """
                                        mkdir -p filtered_contigs
                                        bioawk -c fastx 'length(\$seq) >= 500 { print ">"\$name; print \$seq }' ${contig_file} > filtered_contigs/${x}.fasta
                                        """
}

process metaphylan {
                                        container 'staphb/metaphlan:latest '
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        memory '24 GB'
                                        
                                        

                                        input:
                                        tuple val(x), path(contig_file)

                                        output:
                                        tuple val(x), path("relative_abundance/${x}.txt")


                                        script:
                                        """
                                        mkdir -p relative_abundance
                                       metaphlan ${contig_file} --input_type fasta -o relative_abundance/${x}.txt --tax_lev  a  -t rel_ab
                                        """
}



// Quality control for the assembled genome using Quast

process quast  {
                                        container 'staphb/quast:latest'
                                        publishDir "${params.outdir}", mode: 'copy'


                                        input:
                                        tuple val(x),path(assembled_file)

                                        output:
                                        path("quast_output/${assembled_file.baseName}*")

                                        script:
                                        """
                                        mkdir quast_output
                                        quast.py ${assembled_file}  -o quast_output/${assembled_file.baseName}
                                        """
}

// Creating an index for the draft MAG 

process bowtie_MAGS {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        tuple val(x), path(assembled_file)

                        output:
                        tuple val(x), path("${assembled_file.baseName}.*")
                
                        script:
                        """
                        bowtie2-build ${assembled_file} ${assembled_file.baseName}
                        """
}

// map the indexed genome against the reads 

process bowtie_align_MAGS {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'
                        memory '16 GB'

                        input:
                        tuple val(x), path(paired_reads_1),  path(paired_reads_2),  path(align_input)
                    
                    
                        output:
                        tuple val(x), path("bowtie_output_mags/${x}.sam")
                        
                        script:
                        """
                        mkdir -p bowtie_output_mags
                        
                 
                        bowtie2 -x ${x} -1 ${paired_reads_1} -2 ${paired_reads_2} -S bowtie_output_mags/${x}.sam                   
                        """
}


// Converting sam file to become Bam file followed by sorting and indexing the associated BAM file

process samtools_MAGS {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files/${sam_file.baseName}.sorted.bam"), path("bam_files/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                     mkdir bam_files 
                    samtools view -O BAM ${sam_file} -o  bam_files/${sam_file.baseName}.bam
                    samtools sort  bam_files/${sam_file.baseName}.bam -o  bam_files/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files/${sam_file.baseName}.sorted.bam
                     """
}

//  Quality control statistics for bam file using samtools stats

process samtools_stats {
                                        container 'staphb/samtools:latest'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        containerOptions = "--user root"
                                        cache true

                                        input:
                                        tuple val(x), path(bam_file), path(bam_file_1)

                                        output:
                                        tuple val(x), path("samtools_stats/${bam_file.baseName}.stats")

                                        script:
                                        """
                                        mkdir -p samtools_stats
                                        samtools flagstat ${bam_file} > samtools_stats/${bam_file.baseName}.stats
                                        """
}

process prodigal {
                                        container 'biocontainers/prodigal:v1-2.6.3-4-deb_cv1 '
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        
                                       

                                        input:
                                        tuple val(x), path(filtered_contigs)

                                        output:
                                        tuple val(x), path("prodigal_output/${x}.faa"), path("prodigal_output/${x}.fna"), path("prodigal_output/${x}.gbk") 

                                        script:
                                        """
                                        mkdir -p prodigal_output
                                        prodigal -i ${filtered_contigs} -a prodigal_output/${x}.faa -d prodigal_output/${x}.fna -o prodigal_output/${x}.gbk -f gbk -p meta
                                        """
}


// Determine AMR genes 
process abricate {
					container 'staphb/abricate:latest'
					publishDir "${params.outdir}" , mode: 'copy'
					

					input:
					tuple val(x), path(protein), path(nucleotide), path(gbk)

					output:
					tuple path("abricate_results/${x}.ncbiamr.csv"), 
					path("abricate_results/${x}.resfinderamr.csv"), 
					path("abricate_results/${x}.cardamr.csv"), 
					path("abricate_results/${x}.vfdbamr.csv"), 
					path("abricate_results/${x}.megaresamr.csv"),
					path("abricate_results/${x}.plasmidfinderamr.csv")
					
					script:
					"""
					mkdir abricate_results/
					abricate ${nucleotide} --db ncbi --csv > abricate_results/${x}.ncbiamr.csv 
					abricate ${nucleotide} --db card --csv > abricate_results/${x}.cardamr.csv 
					abricate ${nucleotide} --db resfinder --csv > abricate_results/${x}.resfinderamr.csv  
					abricate ${nucleotide} --db vfdb --csv > abricate_results/${x}.vfdbamr.csv 
					abricate ${nucleotide} --db megares --csv > abricate_results/${x}.megaresamr.csv 
					abricate ${nucleotide} --db plasmidfinder --csv >abricate_results/${x}.plasmidfinderamr.csv
					"""

}



workflow{
my_reads=Channel.fromFilePairs("$params.reads")
my_reference=Channel.fromPath("$params.reference")
//my_reads.view()
fastqc(my_reads)
fastp_ch=fastp(my_reads)
//fastp_ch.view()
trimmomatic_ch=trimmomatic(fastp_ch)
//trimmomatic_ch.view()
fastqc_trimmed(trimmomatic_ch)
bowtie_ch=bowtie_human(my_reference)
//bowtie_ch.view()
fastp_combined_ch=trimmomatic_ch.combine(bowtie_ch)
//fastp_combined_ch.view()
bowtie_aligned_ch=bowtie_align_human(fastp_combined_ch)
//bowtie_aligned_ch.view()
samtools_ch=samtools(bowtie_aligned_ch)
filter_unmapped_reads_ch=filter_unmapped_reads(samtools_ch)
split_file_sort_ch=split_file_sort(filter_unmapped_reads_ch)
split_file_ch=split_file(split_file_sort_ch)
//split_file_ch.view()
spades_ch=spades(split_file_ch)
//spades_ch.view()
bioawk_ch = spades_ch
       .map { x, spadesOutput -> tuple(x, file("${spadesOutput}/${x}/contigs.fasta")) }
 //       .view()
contigs_filtered_ch=bioawk(bioawk_ch)
//rel_abundances_ch=metaphylan(bioawk_ch)
quast_output_ch=quast(contigs_filtered_ch)
bowtie_metagenomes_ch=bowtie_MAGS(contigs_filtered_ch)
bowtie_metagenomes_reads_ch=split_file_ch.combine(bowtie_metagenomes_ch, by: 0)
//bowtie_metagenomes_reads_ch.view()
bowtie_align_metagenome=bowtie_align_MAGS(bowtie_metagenomes_reads_ch)
samtools_metagenome=samtools_MAGS(bowtie_align_metagenome)
samtools_stats(samtools_metagenome)
prodigal_ch=prodigal(contigs_filtered_ch)
//prodigal_ch.view()
abricate(prodigal_ch)

}

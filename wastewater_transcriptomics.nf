#! /usr/bin/nextflow 

params.reads='/home/bioinfo/wastewater_analysis/raw_data/individual_samples/*_{R1,R2}_001.fastq.gz'
params.outdir="/home/bioinfo/wastewater_analysis/results/"
params.reference='/home/bioinfo/wastewater_analysis/raw_data/GRCh38_latest_genomic.fna.gz'
params.rna ='/home/bioinfo/wastewater_analysis/rna_databases/*.fasta'
params.bacteria='/home/bioinfo/wastewater_analysis/rna_databases/silva-bac-16s-id90.fasta'

println "reads: $params.reads"
println "reference: $params.reference"
println "outdir: $params.outdir"
println "rna:     $params.rna"
println "bacteria:  $params.bacteria"



log.info """\
                        Wastewater metatranscriptomic analysis 
                        ------------------------------------------
                        reads:       "$params.reads"
                        reference:   "$params.reference"
                        outdir:      "$params.outdir"
                        rna:          "$params.rna"
                        bacteria:     "$params.bacteria"
                        ------------------------------------------
                        """
                .stripIndent()


// Quality control using fastqc

process fastqc {
                container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
                publishDir "${params.outdir}", mode: 'copy'

                input:
                tuple val(sample_id), path(my_reads_1)

                output:
                tuple val(sample_id), path("*")


                script:
                """
                mkdir -p quality_control
                fastqc ${my_reads_1[0]} ${my_reads_1[1]} -o quality_control

                """

}

// Removal of adapters using fastp

process fastp {
                    container 'quay.io/biocontainers/fastp:0.23.3--h5f740d0_0 '
                    publishDir "${params.outdir}", mode: 'link'

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
                    conda '/home/bioinfo/anaconda3'
                    publishDir "${params.outdir}", mode: 'link'

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
                    
                    trimmomatic  PE ${my_reads_1} ${my_reads_2}  trimmomatic_output/${my_reads_1.baseName}.paired.fastq trimmomatic_output/${my_reads_1.baseName}.unpaired.fastq \
                    trimmomatic_output/${my_reads_2.baseName}.paired.fastq trimmomatic_output/${my_reads_2.baseName}.unpaired.fastq SLIDINGWINDOW:4:20 MINLEN:30
                    """
}

// Checking the quality of trimmed files 

process fastqc_trimmed{
                        container 'biocontainers/fastqc:v0.11.9_cv8'
                        publishDir "${params.outdir}", mode: 'copy'

                        input:
                        tuple val(sample_id), path(trimmed_reads_1), path(trimmed_reads_2), path(trimmed_reads_3), path(trimmed_reads_4)

                        output:
                        tuple val(sample_id), path("*")

                        script:
                        """
                        mkdir -p quality_control_trimmed
                        fastqc ${trimmed_reads_1}  ${trimmed_reads_3}    -o quality_control_trimmed

                        """
}

// merge forward and reverse reads using PEAR 



process PEAR {
                  container 'quay.io/biocontainers/pear:0.9.6--h67092d7_9 '
                  publishDir "${params.outdir}", mode: 'copy'

            
                  input:
                  tuple val(x), path(paired_reads_1), path(unpaired_reads_1), path(paired_reads_2), path(unpaired_reads_2)

                  output:
                  tuple val(x), path("*")

                  script:
                  """
                  mkdir -p merged_reads
                  pear -f ${paired_reads_1} -r ${paired_reads_2} -o merged_reads/${x}
                  """
}




// Creating an index for the reference genome 

process bwa_index {
                        conda '/home/bioinfo/anaconda3'
                        publishDir "${params.outdir}", mode: 'link'
                        memory '28 GB'

                        input:
                        path reference

                        output:
                        tuple val("GRCh38_latest_genomic.fna.gz"), path("*")
                
                        script:
                        """
                        bwa index ${reference} 
                        """
}

// map the indexed genome against the reads 

process bwa_align_human {
                        conda '/home/bioinfo/anaconda3'
                        publishDir "${params.outdir}", mode: 'link'
                        memory '28 GB'

                        input:
                        tuple val(index), path(indexes), val(x),  path(reads_1), path(reads_2), path(reads_3), path(reads_4)
                    
                    
                        output:
                        tuple val(x), path("bam_files/${x}.bam")
                        
                        script:
                        """
                        mkdir -p bam_files indexed_files 
                        mv ${indexes} indexed_files/
                                        
                        bwa mem indexed_files/${index}  ${reads_1} ${reads_3} | samtools view -bS - > bam_files/${x}.bam  
                                       
                        """
}


// Converting sam file to become Bam file followed by sorting and indexing the associated BAM file

process samtools {
                    conda '/home/bioinfo/anaconda3'
                     publishDir "${params.outdir}" , mode: 'link'
                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files_sorted/${sam_file.baseName}.sorted.bam"), path("bam_files_sorted/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                    mkdir -p bam_files_sorted 
                    samtools sort ${sam_file} -o  bam_files_sorted/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files_sorted/${sam_file.baseName}.sorted.bam
                    
                     """
}

// filter out unmapped reads 
process filter_unmapped_reads {
                    conda '/home/bioinfo/anaconda3'
                     publishDir "${params.outdir}" , mode: 'link'
                    

                     input:
                     tuple val(x),path(bam_file), path(bam_file_index)

                     output:
                     tuple val(x), path("bam_files_unmapped/${bam_file.baseName}.unmapped.bam"), path("bam_files_unmapped/${bam_file.baseName}.unmapped.sorted.bam")

                     script:
                     """
                     mkdir -p bam_files_unmapped
                     samtools view -b -f 4 ${bam_file} -o bam_files_unmapped/${bam_file.baseName}.unmapped.bam
                     samtools sort bam_files_unmapped/${bam_file.baseName}.unmapped.bam -o bam_files_unmapped/${bam_file.baseName}.unmapped.sorted.bam -O BAM
                    
                     """
}


// split files to fastq 

process split_file  {
                     conda '/home/bioinfo/anaconda3'
                     publishDir "${params.outdir}" , mode: 'link'
                    

                     input:
                     tuple val(x),path(bam_file_1), path(bam_file_2)

                     output:
                     tuple val(x), path("bam_files_splitted/${x}_1.fastq"), path("bam_files_splitted/${x}_2.fastq")

                     script:
                     """
                     mkdir -p bam_files_splitted
                     samtools fastq ${bam_file_2} -1  bam_files_splitted/${x}_1.fastq -2 bam_files_splitted/${x}_2.fastq
                           
                     """
}


process SORTMERNA {

                  container 'quay.io/biocontainers/sortmerna:4.3.6--h9ee0642_0'
                  publishDir "${params.outdir}", mode: 'link'
                  memory '24 GB'

            
                  input:
                  tuple val(x), path(reads_1), path(reads_2), path(database_1), path(database_2), path(database_3), path(database_4), path(database_5), path(database_6)

                  output:
                  tuple val(x), path("*")

                  script:
                  """
                  mkdir -p sortmeresults 
                  sortmerna  --ref ${database_1} --ref ${database_2} --ref ${database_3} --ref ${database_4} --ref ${database_5} --ref ${database_6}  --reads ${reads_1} --reads ${reads_2} --workdir sortmeresults --fastx --aligned --other
                  """
}

// assmble transcriptome using trinity 

process trinity {
                                        container 'quay.io/biocontainers/trinity:2.15.1--pl5321h146fbdb_3'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        memory '24 GB'
                                        

                                        input:
                                        tuple val(x), path(trimmed_reads_1)

                                        output:
                                        tuple val(x), path("*")


                                        script:
                                        """
                                        mkdir -p trinity_out_dir 
                                        Trinity  --seqType fq  --single ${trimmed_reads_1} --output trinity_out_dir/trinity_${x} --max_memory 24G  --full_cleanup 
                                        """
}

// Running de novo assembly for the samples using   SPADES

process spades {
                                        container 'quay.io/biocontainers/spades:3.15.5--h95f258a_0  '
                                        publishDir "${params.outdir}" , mode: 'link'
                                        memory '28 GB'
                                        

                                        input:
                                        tuple val(x), path(trimmed_reads_1)

                                        output:
                                        tuple val(x), path("*")


                                        script:
                                        """
                                        mkdir spades_output
                                        spades.py -s ${trimmed_reads_1}  --meta  -o spades_output/${x} --tmp-dir . --threads 1
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
                        tuple val(x), path("${x}.*")
                
                        script:
                        """
                        bowtie2-build ${assembled_file} ${x}
                        """
}

// map the indexed genome against the reads 

process bowtie_align_MAGS {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'
                        memory '16 GB'

                        input:
                        tuple val(x), path(paired_reads_1),  path(align_input)
                    
                    
                        output:
                        tuple val(x), path("bowtie_output_mags/${x}.sam")
                        
                        script:
                        """
                        mkdir -p bowtie_output_mags
                        
                 
                        bowtie2 -x ${x} -U ${paired_reads_1} -S bowtie_output_mags/${x}.sam                   
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

                                        input:/trinity_out_dir/
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


// Creating an index for the reference genome 

process bowtie_bacteria {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'
                        memory '16 GB'

                        input:
                        path reference

                        output:
                        tuple val("bowtie_bacteria"), path("bowtie_bacteria*")
                
                        script:
                        """
                        bowtie2-build ${reference}  bowtie_bacteria
                        """
}

// map the indexed genome against the reads 

process bowtie_align_bacteria {
                        container 'biocontainers/bowtie2:v2.4.1_cv1'
                        publishDir "${params.outdir}", mode: 'copy'
                        memory '16 GB'

                        input:
                        tuple val(x), path(paired_reads_1), val(index), path(align_input)
                    
                    
                        output:
                        tuple val(x), path("bowtie_output_bacteria/${x}.sam")
                        
                        script:
                        """
                        mkdir -p bowtie_output_bacteria
                                        
                        bowtie2 -x ${index} -U ${paired_reads_1} -S bowtie_output_bacteria/${x}.sam                   
                        """
}


// Converting sam file to become Bam file followed by sorting and indexing the associated BAM file

process samtools_bacteria {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(sam_file)

                     output:
                     tuple val(x), path("bam_files_bacteria/${sam_file.baseName}.sorted.bam"), path("bam_files_bacteria/${sam_file.baseName}.sorted.bam.bai")

                     script:
                     """
                     mkdir bam_files_bacteria 
                    samtools view -O BAM ${sam_file} -o  bam_files_bacteria/${sam_file.baseName}.bam
                    samtools sort  bam_files_bacteria/${sam_file.baseName}.bam -o  bam_files_bacteria/${sam_file.baseName}.sorted.bam -O BAM
                    samtools index  bam_files_bacteria/${sam_file.baseName}.sorted.bam
                     """
}


process samtools_stats_bacteria {
                                        container 'staphb/samtools:latest'
                                        publishDir "${params.outdir}" , mode: 'copy'
                                        containerOptions = "--user root"
                                        cache true

                                        input:
                                        tuple val(x), path(bam_file), path(bam_file_1)

                                        output:
                                        tuple val(x), path("samtools_stats_bacteria/${bam_file.baseName}.stats")

                                        script:
                                        """
                                        mkdir -p samtools_stats_bacteria
                                        samtools flagstat ${bam_file} > samtools_stats_bacteria/${bam_file.baseName}.stats
                                        """
}


process mapped_bacteria {
                     container 'staphb/samtools:latest'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x),path(bam_file), path(bai_file)

                     output:
                     tuple val(x), path("bacteria_mapped/${x}.fasta")

                     script:
                     """
                     mkdir -p bacteria_mapped
                     samtools fasta -F 4 ${bam_file} > bacteria_mapped/${x}.fasta
                           
                     """
}


process vsearch {
                     container 'biocontainers/vsearch:v2.10.4-1-deb_cv1 '
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x), path(fasta_file)

                     output:
                     tuple val(x), path("dereplicate_results/${x}.fasta")

                     script:
                     """
                     mkdir -p dereplicate_results
                     vsearch --derep_fulllength ${fasta_file} -sizeout -output dereplicate_results/${x}.fasta
                           
                     """
}


process chimera_detection {
                     container 'biocontainers/vsearch:v2.10.4-1-deb_cv1'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x), path(fasta_file)

                     output:
                     tuple val(x), path("non_chimera/${x}.fasta")

                     script:
                     """
                     mkdir -p non_chimera
                     vsearch  --uchime_denovo  ${fasta_file}  --abskew   1.5   --nonchimeras   non_chimera/${x}.fasta --fasta_width 0
                           
                     """
}


process greedy_delimitation {
                     container 'biocontainers/vsearch:v2.10.4-1-deb_cv1'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x), path(fasta_file)

                     output:
                     tuple val(x), path("greedy_delimit/${x}.otu.fasta"), path("greedy_delimit/${x}.uc")

                     script:
                     """
                     mkdir -p greedy_delimit
                     vsearch  --cluster_size  ${fasta_file} --id 0.97 --centroids   greedy_delimit/${x}.otu.fasta --sizein --relabel otu --uc greedy_delimit/${x}.uc
                           
                     """
}

process mapped_reads_otus {
                     container 'biocontainers/vsearch:v2.10.4-1-deb_cv1'
                     publishDir "${params.outdir}" , mode: 'copy'
                    

                     input:
                     tuple val(x), path(merged_file), path(otu_file), path(uc_file)
                     

                     output:
                     tuple val(x), path("otu_counts/${x}.tsv")

                     script:
                     """
                     mkdir -p otu_counts
                     vsearch --usearch_global ${merged_file} -db ${otu_file} -id 0.97 -otutabout otu_counts/${x}.tsv
                           
                     """
}


workflow{

/// Read the different co-assembled reads and databases 
my_reads=Channel.fromFilePairs("$params.reads")
my_reference=Channel.fromPath("$params.reference")
my_rna_databases=Channel.fromPath("$params.rna")
bacteria=Channel.fromPath("$params.bacteria")
//my_reads.view()
//my_rna_databases.view()
rna_databases=my_rna_databases.collect()
//rna_databases.view()


/// Check the quality of the sequences using FASTQC 
fastqc(my_reads)

/// Remove adapters using fastP 
fastp_ch=fastp(my_reads)
//fastp_ch.view()

/// Trim low quality sequences using Trimmomatic 
trimmomatic_ch=trimmomatic(fastp_ch)
//trimmomatic_ch.view()

//// Check the quality of trimmed reads 
fastqc_trimmed(trimmomatic_ch)


//// Merge the trimmed reads 
//merged_fastq_files=PEAR(trimmomatic_ch)

/// select only the assembled merged reads 
//merged_ch=merged_fastq_files.map{x, merged_reads -> tuple(x, file("${merged_reads}/*.assembled.fastq")) }
//merged_ch.view()

/// index the human reference genome 
bowtie_ch=bwa_index(my_reference)
//bowtie_ch.view()


/// combine merged reads with indexed reads 
bowtie_combined_ch=bowtie_ch.combine(trimmomatic_ch)
//bowtie_combined_ch.view()

/// alig reads with the human genome 
bowtie_aligned_ch=bwa_align_human(bowtie_combined_ch)
//bowtie_aligned_ch.view()

//// convert sam file from mapping to become  bam file 
samtools_ch=samtools(bowtie_aligned_ch)
//samtools_ch.view()


//// Filter human reads to remain with non-host reads 
filter_unmapped_reads_ch=filter_unmapped_reads(samtools_ch)
//filter_unmapped_reads_ch.view()

/// split the file and convert bam file to fastq file 
//split_file_sort_ch=split_file_sort(filter_unmapped_reads_ch)
//split_file_sort_ch.view()
split_file_ch=split_file(filter_unmapped_reads_ch)
//split_file_ch.view()

//// combine the reads with the baxterial silva databases 
rna_remove=split_file_ch.combine(rna_databases)
//rna_remove.view()

/// use sortmeRNA to remove rRNA reads from the sample 
non_rna_reads=SORTMERNA(rna_remove)
non_rna_reads.view()

/// select other reads to use in downstrean processess 
cleaned_reads = non_rna_reads.map { x , cleaned_reads -> tuple(x, file("${cleaned_reads}/out/other.fq")) }
cleaned_reads.view()

//// Assemble reads using trinity 

trinity=trinity(cleaned_reads)
//trinity.view()
//trinity_output =trinity.map{ x, trinity_output -> tuple (x, file("${trinity_output}/*.fasta"))}
//trinity_output.view()

///// assess the quality of the reads 
//quast_output_ch=quast(trinity_output)

/// de novo assembly using spades 
spades_ch=spades(cleaned_reads)
//spades_ch.view()
//bioawk_ch = spades_ch
  //     .map { x, spadesOutput -> tuple(x, file("${spadesOutput}/${x}/transcripts.fasta")) }
   //   .view()
//contigs_filtered_ch=bioawk(bioawk_ch)
//rel_abundances_ch=metaphylan(bioawk_ch)
//quast_output_ch=quast(contigs_filtered_ch)

// index trinity output 
//bowtie_metagenomes_ch=bowtie_MAGS(trinity_output)


//bowtie_metagenomes_reads_ch=cleaned_reads.combine(bowtie_metagenomes_ch, by: 0)
//bowtie_metagenomes_reads_ch.view()
// align reads to indexed genome 
//bowtie_align_metagenome=bowtie_align_MAGS(bowtie_metagenomes_reads_ch)

// convert sam file to bam file 
//samtools_metagenome=samtools_MAGS(bowtie_align_metagenome)

// check mapped statistics 
//samtools_stats(samtools_metagenome)

// use prodigal to predict open reading frame 
//prodigal_ch=prodigal(trinity_output)
//prodigal_ch.view()

// use abricate to predict resistance 
//abricate(prodigal_ch)





// Analysis of the files which aligned to the silva database. Specifically mapping against the bacterial-16s rRNA to classify abundance of different reads 

/// getting the aligned sequences from sortmeRNA 
//aligned_rna= non_rna_reads.map { x , aligned_rna -> tuple(x, file("${aligned_rna}/out/aligned.fq")) }
//aligned_rna.view()

/// Index the bacterial genome 16s rRNA database 
//bacteria_index=bowtie_bacteria(bacteria)
//bacteria_index.view()

/// combine each of the bacteria aligned files with the bacterial indexed genome 
//bacteria_combined=aligned_rna.combine(bacteria_index)
//bacteria_combined.view()

// Align the reads against the indexed bacterial geneome 
//bacteria_aligned=bowtie_align_bacteria(bacteria_combined)

// convert the sam file to become bam files 
//bam_file_bacteria=samtools_bacteria(bacteria_aligned)
//bam_file_bacteria.view()


///  Check the percentages of mapped reads 
//samtools_stats_bacteria(bam_file_bacteria)

/// filter out the mapped bacterial reads 
//bacteria_mapped_ch=mapped_bacteria(bam_file_bacteria)

/// dereplicate mapped reads to remains with non redudant reads 
//dereplication=vsearch(bacteria_mapped_ch)

/// remove chimeras 
//chimera_ch=chimera_detection(dereplication)

///  use greedy delimiatation to identify OTUs 
//greedy_delimit_ch=greedy_delimitation(chimera_ch)

/// Get the reads which mount to the otus 
//merged_greedy=bacteria_mapped_ch.combine(greedy_delimit_ch, by:0)

//// get the counts table
//merged_greedy.view()
//otus_mapped_ch=mapped_reads_otus(merged_greedy)
}


#! /usr/bin/nextflow 

params.reads='/home/bioinfo/wastewater_analysis/raw_data/individual_samples/*_{R1,R2}_001.fastq.gz'
params.outdir="/home/bioinfo/wastewater_analysis/results/"
params.reference='/home/bioinfo/wastewater_analysis/raw_data/GRCh38_latest_genomic.fna.gz'
params.rna ='/home/bioinfo/wastewater_analysis/rna_databases/*.fasta'

println "reads: $params.reads"
println "reference: $params.reference"
println "outdir: $params.outdir"
println "rna:     $params.rna"




log.info """\
                        Wastewater metatranscriptomic analysis 
                        ------------------------------------------
                        reads:       "$params.reads"
                        reference:   "$params.reference"
                        outdir:      "$params.outdir"
                        rna:          "$params.rna"
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

}

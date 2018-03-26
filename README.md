# RunGBS
bash wrapper for running GBS/ddRAD pipeline

**Download links for the tools:**
1. bowtie2 : https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1
2. Parallel : https://www.gnu.org/software/parallel/
3. sabre : https://github.com/najoshi/sabre
4. cutadapt : Recommended to install through bioconda (http://cutadapt.readthedocs.io/en/stable/installation.html)
5. bcftools : Recommended to install through bioconda (http://www.htslib.org/download/)
6. samtools : Recommended to install through bioconda (http://www.htslib.org/download/)



```
        USAGE: runGBS.pl --se/--pe -r <reference.fa> -k <GBS_KeyFile.txt> -i <fastq_dir_path>                                                                                                
                                                                                                                                                                                             
        Options:                                                                                                                                                                             
        --se|S  - Single-end data                                                                                                                                                            
        --pe|P  - Paired-end data                                                                                                                                                            
        --reference|-r  - Reference sequence file. (.fasta or .fa)                                                                                                                           
        --keyfile|-k    - Path for Keyfile.                                                                                                                                                  
        --fastq|-i      - Path for fastq directory (fastq files should named as FLOWCELLIF_LANE.fq.gz)                                                                                       
                                                                                                                                                                                             
        Fastq processing:                                                                                                                                                                    
        --adapter|-a    - Adapter sequence to use for cutadapt tool. (default: GAGATCGGAAGAGCGGG)                                                                                            
        --minlength|-m  - Minimum read length to be considered. (default: 50)

        Performance:
        --threads|-p    - Number of threads/cores to be used for bowtie2/samtools/bcftools (default: 2)
        --npcore|-n     - Number of threads/cores to be used for parallel (default: 2)

        Tools:
        --sabre - Paht for sabre executable (default: $PATH/sabre)
        --cutadapt      - Paht for cutadapt executable (default: $PATH/cutadapt)
        --bowtie2       - Paht for bowtie2 executable (default: $PATH/bowtie2)
        --samtools      - Paht for samtools executable (default: $PATH/samtools)
        --bcftools      - Paht for bcftools executable (default: $PATH/bcftools)

```

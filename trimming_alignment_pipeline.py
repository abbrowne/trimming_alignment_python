## ----------------================Directions for using this script================----------------
## Ensure that fastq.gz files are in one directory and pass this script the following arguments:
##  1: Working directory   2: Genome directory for alignment   3: GTF file path for alignment
## * You will need to edit the parsing of the gz files in the script or edit the names of the *
##       * files themselves to match the following format: sample_name_R#.fastq.gz *
##                 * with # being 1 and 2 for each read pair member *
##      * e.g. sample_A_R1.fastq.gz and sample_A_R2.fastq.gz would be one read pair *

import os
import re
import pandas as pd
import numpy as np
import sys


##Assign main directory for trimming and alignment
main_dir = sys.argv[1]
genome_dir = sys.argv[2]
gtf_dir = sys.argv[3]

os.chdir(main_dir)
##Get gz_files
gz_files = pd.Series(os.listdir('.'))
gz_files = gz_files[gz_files.str.contains('.fastq.gz')].sort_values()
gz_files = gz_files.reset_index(drop=True)

##Strip ending from file names and parse into base names of reads from each pair and main names of the combined pair
temp_base_names = pd.Series([])
temp_main_names = pd.Series([])
for base_name in gz_files:
        temp_base_name = pd.Series([base_name.replace('.fastq.gz','')])
        temp_base_names = temp_base_names.append(temp_base_name)
        temp_main_name = pd.Series([base_name.replace('_R1.fastq.gz','')])
        temp_main_names = temp_main_names.append(temp_main_name)

base_names = temp_base_names.copy().reset_index(drop=True)
main_names = temp_main_names[::2].copy().sort_values().reset_index(drop=True)

##Check file names
gz_files.head()
base_names.head()
main_names.head()

##Initialize text and counters for full trimming run
full_trimming_run = '#!/bin/bash\n'
R1_temp = 0
R2_temp = 1

##Run trimming pipeline on each combined read pair
for main_name in main_names:
        full_trimming_run = full_trimming_run + 'bsub < ./' + main_name + '_trimming.sh\nsleep 2\n'
        ##Initialize trimming subrun text
        trimming_run = '''#!/bin/bash
        #BSUB -J {}_trim
        #BSUB -n 1
        #BSUB -R "rusage[mem=8000]"
        #BSUB -R "span[ptile=1]"
        #BSUB -W 04:00
        #BSUB -o %J.stdout
        #BSUB -eo %J.stderr
        #BSUB -L /bin/bash

        echo Started run at:
        date
        cd {}
        module load java
        module load fastqc
        module load trim_galore
        fastqc -o {} {}
        fastqc -o {} {}
        echo Finished pretrim QC at:
        date
        trim_galore --paired --illumina {} {}
        echo Finished trimming at:
        date
        fastqc -o {} {}_val_1.fq.gz
        fastqc -o {} {}_val_2.fq.gz
        echo Finished posttrim QC at:
        date
        echo Finished run at:
        date
        '''.format(main_name, main_dir, main_dir, gz_files[R1_temp], main_dir, gz_files[R2_temp], gz_files[R1_temp], gz_files[R2_temp], main_dir, base_names[R1_temp], main_dir, base_names[R2_temp])
        R1_temp += 2
        R2_temp += 2
        ##Write trimming subrun job script to file
        open_file = open(str(main_name + '_trimming.sh'), 'w')
        open_file.write(trimming_run)
        open_file.close()

##Write full trimming run job submission to file
full_trim = open(str('full_trimming_run.sh'), 'w')
full_trim.write(full_trimming_run)
full_trim.close()
##Provide full user access to sh files
os.system('chmod u+wrx *.sh')

##Initialize text and counters for full trimming run
full_alignment_run = '#!/bin/bash\n'
R1_temp = 0
R2_temp = 1

##Run alignment pipeline on each combined trimmed read pair
for main_name in main_names:
        full_alignment_run = full_alignment_run + 'bsub < ./' + main_name + '_alignment.sh\nsleep 2\n'
        ##Initialize trimming subrun text
        alignment_run = '''#!/bin/bash
        #BSUB -J {}_align
        #BSUB -n 4
        #BSUB -R "rusage[mem=8000]"
        #BSUB -R "span[hosts=1]"
        #BSUB -W 04:00
        #BSUB -o %J.stdout
        #BSUB -eo %J.stderr
        #BSUB -L /bin/bash

        echo Started run at:
        date
        cd {}
        module load star/2.4.0j
        module load samtools
        module load java
        module load picard/1.112
        module load pigz
        module load subread

        ##Perform alignment with STAR
        STAR --genomeDir {} --sjdbGTFfile {} --runThreadN 4 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMmode Full --outFileNamePrefix {}{}. --readFilesCommand zcat --readFilesIn {}_val_1.fq.gz,{}_val_2.fq.gz
        echo Finished alignment at:
        date
        mv {}{}.Aligned.sortedByCoord.out.bam {}.accepted_hits.sort.coord.bam
        ##Build bam index
        java -jar $PICARD_HOME/BuildBamIndex.jar I={}{}.accepted_hits.sort.coord.bam
        cp {}{}.accepted_hits.sort.coord.bai {}{}.accepted_hits.sort.coord.bam.bai
        ##Run featureCounts on bam files
        featureCounts -T 10 -t exon -g transcript_id -p -a {} -o {}.exon.transcriptID.txt {}{}.accepted_hits.sort.coord.bam
        featureCounts -T 10 -t exon -g gene_id -p -a {} -o {}.exon.geneID.txt {}{}.accepted_hits.sort.coord.bam
        featureCounts -T 10 -t exon -f -O -p -a {} -o {}.exon.txt {}{}.accepted_hits.sort.coord.bam
        featureCounts -T 10 -t exon -g gene_id --primary -O -p -a {} -o {}.primary.txt {}{}.accepted_hits.sort.coord.bam
        echo Finished featureCounts at:
        date
        '''.format(main_name, main_dir, genome_dir, gtf_dir, main_dir, main_name, base_names[R1_temp], base_names[R2_temp], main_dir, main_name, main_name, main_dir, main_name, main_dir, main_name, main_dir, main_name, gtf_dir, main_name, main_dir, main_name, gtf_dir, main_name, main_dir, main_name, gtf_dir, main_name, main_dir, main_name, gtf_dir, main_name, main_dir, main_name)
        R1_temp += 2
        R2_temp += 2
        ##Write alignment subrun job script to file
        open_file = open(str(main_name + '_alignment.sh'), 'w')
        open_file.write(alignment_run)
        open_file.close()

##Write full alignment run job submission to file
full_align = open(str('full_alignment_run.sh'), 'w')
full_align.write(full_alignment_run)
full_align.close()
os.system('chmod u+wrx *.sh')


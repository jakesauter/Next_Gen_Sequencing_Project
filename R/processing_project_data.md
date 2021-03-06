---
title: "First-Steps: Processing Project Data"
author: "Jake Sauter"
date: "2/17/2021"
output: 
  html_document: 
    toc: true
    keep_md: true
editor_options: 
  chunk_output_type: console
---




## **Context of Project FastQ Files**

Information the `fastq` files used in this study can be found [here](https://www.ncbi.nlm.nih.gov/Traces/study/?page=3&acc=PRJNA515044&o=diagnosis_sam_s%3Ad%253Bacc_s%3Bacc_s%3Aa) 
as well as [NCBI GEO](https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE125050)

**Where did you get it from?**: NCBI SRA Run Selector.

**What publication is it linked to?**: 
[Alzheimer's Patient Microglia Exhibit Enhanced Aging and Unique Transcriptional Activation](https://pubmed.ncbi.nlm.nih.gov/32610143/)

**Who generated the data?**: Genentech, Inc.

**How was the RNA extracted?**: 

> <font size=2.5> For whole tissue RNA studies (GSE125583), frozen tissue was sectioned in approximately 8 slices 40mm thick and stored at 80 degrees C. Tissue was homogenized in 1 mL QIAzol with 5 mm stainless steel beads using a Tissuelyzer (20 Hz for 4 min). After homogenization,200mL of choloroform were added to the cleared lysate (1 min at 12,000 rcf.at 4 degrees C), vigorously shook and incubated at room tem-perature 2-3 min. Samples were centrifuged for 15 min at 12,000 rcf.at 4 degrees C and the upper aqueous phase was transferred to a newtube. **RNA was extracted using QIAGEN miRNeasy mini columns**, yielding samples with RNA integrity (RIN) scores averaging 6.5. </font>

![](/home/x1/Documents/Weill_Cornell/ANGSD/angsd_project/imgs/rna_extraction_step.png)

**What library prep was used?**: 

> <font size=2.5> Given the highly fragmented condition of our sorted cell RNA preps, we chose the NuGEN Ovation RNA-Seq System V2 kit for cDNA synthesis since it uses random oligos for cDNA priming. We knew this would result in high percentages of intronic and non-coding RNA reads, but our priority was to sample across all exons instead of having an extreme 3' bias and reduced complexity in our library. (Only exonic reads were counted toward nRPKM values.) Generated cDNA was sheared to 150-200bp size using LE220 ultrasonicator (Covaris). Following shearing, the size of cDNA was determined by Bioanalyzer DNA 1000 Kit (Agilent) and quantity was determined by Qubit dsDNA BR Assay (Life Technologies). Sheared cDNA was subjected to **library generation, starting at end repair step, using Illumina’s TruSeq RNA Sample Preparation Kit v2 (Illumina)**. Size of the libraries was confirmed using 4200 TapeStationand High Sensitivity D1K screen tape (Agilent Technologies) and their concentration was determined using KAPA Library Quantifi-cation kits. **The libraries were multiplexed within cell types and then sequenced on Illumina HiSeq2500 (Illumina) to generate 50M of single end 50bp reads**. </font>

**What cell type was used?**: The cell types used in this data include Astrocytes, Endothelial cells, Microglia and Neurons.

**What was the treatment/experimental condition?**: Tissues originating from either healthy control condition brains or brains of diseased patients known to possess Alzheimer's Disease. 

**What sequencing platform was used?**: Illumina HiSeq 2500
    
## **Downloading the Fastq Files**

In order to first view the run accession metadata file for 
the fastq files of interest: 

  * Navigate to: 
      https://www.ncbi.nlm.nih.gov/Traces/study/?page=3&acc=PRJNA515044&o=diagnosis_sam_s%3Ad%3Bacc_s%3Aa
  
  * Download "metadata" file
  * Run the code below to download

Since I have not found a way to do this in an automated fashion, I have
downloaded the file and am hosting it myself through a Github repository
made to host this assignment. Thus the **above steps do not need to be 
followed and the command shown below can just be run on the target system**.


```bash
mkdir accession_tables

wget https://raw.githubusercontent.com/jakesauter/Next_Gen_Sequencing_Project/main/accession_tables/Sra_run_table.txt --output-document="accession_tables/Sra_run_table.txt"
```

**Exploring the dataset**

Lets explore the distribution of experimental conditions in this dataset.


```r
library(dplyr)
library(knitr)
```



```r
filepath <- "accession_tables/Sra_run_table.txt"

df <- read.csv(filepath, 
               header = TRUE, 
               stringsAsFactors = FALSE) %>% 
  as.data.frame()

df %>% 
  group_by(Diagnosis, Cell_type) %>% 
  summarise(Num_Samples=n()) %>% 
  kable()
```



|Diagnosis |Cell_type   | Num_Samples|
|:---------|:-----------|-----------:|
|AD        |astrocyte   |           7|
|AD        |endothelial |          10|
|AD        |myeloid     |          10|
|AD        |neuron      |          21|
|Control   |astrocyte   |          12|
|Control   |endothelial |          17|
|Control   |myeloid     |          15|
|Control   |neuron      |          21|

We can see that the 4 main cell types are represented across
control and Alzheimer's Disease (AD). Lets choose only myeloid cells moving forward, and explore the difference in sex of collected tissues.


```r
df %>% 
  filter(Cell_type == 'myeloid', 
         !is.na(sex)) %>% 
  group_by(Diagnosis, sex) %>% 
  summarise(Num_Samples=n()) %>% 
  kable()
```



|Diagnosis |sex    | Num_Samples|
|:---------|:------|-----------:|
|AD        |female |           3|
|AD        |male   |           5|
|Control   |female |           4|
|Control   |male   |           9|

From above we see that our smallest category for myeloid cells is AD female with a group size of 3. Moving forward, we will strive to collect 3 samples from every group for initial quality control and mapping testing.


```r
# Initially, these accession numbers were randomly sampled
# from the df data frame below with sample_n(3), though string
# literals are used now in order for this document to produce
# the same runs upon rendering
acc_nums <- c('SRR8440514', 'SRR8440538', 'SRR8440448', 'SRR8440501', 'SRR8440463', 'SRR8440477', 'SRR8440443', 'SRR8440511', 'SRR8440518', 'SRR8440482', 'SRR8440524', 'SRR8440539')


df <- df %>% 
  filter(Cell_type == 'myeloid', 
         !is.na(sex), 
         Run %in% acc_nums) %>% 
  group_by(Diagnosis, sex) %>%
  ungroup() %>% 
  select("Run", "Cell_type", "sex",  "Diagnosis", "expired_age")

df %>% 
  kable()
```



|Run        |Cell_type |sex    |Diagnosis |expired_age |
|:----------|:---------|:------|:---------|:-----------|
|SRR8440477 |myeloid   |male   |AD        |88          |
|SRR8440501 |myeloid   |male   |AD        |79          |
|SRR8440511 |myeloid   |female |Control   |59          |
|SRR8440518 |myeloid   |female |Control   |88          |
|SRR8440443 |myeloid   |female |Control   |>90         |
|SRR8440448 |myeloid   |female |AD        |84          |
|SRR8440482 |myeloid   |male   |Control   |>90         |
|SRR8440514 |myeloid   |female |AD        |78          |
|SRR8440524 |myeloid   |male   |Control   |86          |
|SRR8440538 |myeloid   |female |AD        |84          |
|SRR8440539 |myeloid   |male   |Control   |>90         |
|SRR8440463 |myeloid   |male   |AD        |72          |

We can then use this accession list to download our 
files of interest.


```bash
cell_type="myeloid_cells"
acc_nums=(SRR8440514 SRR8440538 SRR8440448 SRR8440501 SRR8440463 SRR8440477 SRR8440443 SRR8440511 SRR8440518 SRR8440482 SRR8440524 SRR8440539)

for acc_num in "${acc_nums[@]}" ; do 
  fasterq-dump --outdir $cell_type $acc_num
  gzip ${cell_type}/${acc_num}.fastq
done
```

## **FastQC Quality Control**


```bash
cell_type="myeloid_cells"
fastq_files=$(ls "${cell_type}"/*.fastq.gz)
fastqc_outdir="${cell_type}_fastqc"

mkdir "${cell_type}_fastqc"

for fastq_file in $fastq_files; do
  fastqc $fastq_file --noextract \
    --outdir=$fastqc_outdir \
    --threads=12
done
```



```bash
cd $fastqc_outdir
multiqc . 
```

The above command will generate a `multiqc_report.html` file which 
we can inspect to determine the overall quality of the sequenced 
libraries. 

### **Per Base Sequence Quality**

Overall sequence quality of the sequenced libraries appears to be very 
high when view the **Per Base Sequence Quality** score in the resultant
`multiqc_report.html`. 

![](/home/x1/Documents/Weill_Cornell/ANGSD/angsd_project/imgs/fastqc_per_base_sequence_quality_plot.png){}

### **Over-Represented Sequences**

`MultiQC` output has also informed us about some warnings
regarding over-represented sequences. Unfortunality, this
sequence is not shown in `MultiQC` output, so we will 
inspect a `FastQC` output html file to gain more information 
about this. For this exercise I have chosen to view `SRR8440524_fastqc.html`, 
which shows us below that we may have some Illumina TruSeq Adapter 
dimerization contamination. 

![](/home/x1/Documents/Weill_Cornell/ANGSD/angsd_project/imgs/overp_seqs.png)

### **TrimGalore**

Due to the highly frequent over-represented sequences seen above, I have opted
to try to use [`TrimGalore`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) on these reads to remove the potential adapter sequence.

Notice below that I did not choose to use the `--fastqc` or `--fastqc_args`
parameters to `trim_galore`. I chose not to use these options as 
in order to parallelize `FastQC`, it must be working on multiple files
at a time. Thus, in order to more efficiently make use of the system it would
actually be quicker to run `FastQC` independent of `trim_galore` afterwards.


**FASTQC On Trimmed Sequences**


```bash
cell_type="myeloid_cells"
fastq_files=$(ls "${cell_type}"/*.fastq.gz)

trim_galore_outdir="trimmed_${cell_type}"
fastqc_outdir="trimmed_${cell_type}_fastqc"
mkdir $fastqc_outdir

trim_galore $fastq_files \
  --output_dir="${trim_galore_outdir}" \
  --cores 8 --phred33 --gzip
  
trimmed_files=$(ls $trim_galore_outdir/*.fq.gz)
fastqc $trimmed_files \
  --outdir=$fastqc_outdir \
  --threads=12 --noextract
  
cd $fastqc_outdir
multiqc . 
```

**Trim Reports**

Below we can see that `TrimGalore` did not remove too many 
sequences from any one file, with the largest amount of sequences
being **2.3%**, coming from sample **SRR8440443**

![](processing_project_data_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

**Code for above plot**


```bash
tail -n 2 $trim_galore_outdir/*_trimming_report.txt > trimming_reports.txt
```


```r
library(stringr)
library(magrittr)
library(ggplot2)

lines <- readLines('trimming_reports.txt')
lines <- lines[lines != ""]

removed_seq_data <- vector('numeric', length(lines)/2)

for (i in seq(1,length(lines), 2)) {
  sample_name <- lines[i] %>% 
    str_extract("SRR[0-9]*")
  removed_seq_percentage <- lines[i+1] %>% 
    str_extract("[0-9]\\.[0-9]")
  removed_seq_data[(i+1)/2] <- removed_seq_percentage
  names(removed_seq_data)[(i+1)/2] <- sample_name
}

removed_seq_data %>% 
  data.frame(Sample=names(.), 
             Percent_Dropped=., 
             row.names = NULL) %>% 
  ggplot(aes(x=Sample, y=Percent_Dropped)) +
    geom_col() + 
    ylab("% Sequenced Dropped (TrimGalore)") + 
    xlab("") + 
    theme(axis.text.x=element_text(angle = -25, face='bold'), 
          axis.text.y=element_text(face='bold'))
```


As well as total sequences removed, we can see that the total sequence
length distribution has not been devastated after the use of `TrimGalore`, 
which should give us some trust that it has not removed large chunks 
from too many sequences.

![](/home/x1/Documents/Weill_Cornell/ANGSD/angsd_project/imgs/post_trim_seq_len_dist.png){width=50%}


### **Per Base Sequence Content**

Some issues were spotted in the `FastQC` reports in the per base
sequence content section for some of the samples. A few of the samples
seem to show a cyclic pattern on this plot, indicating a potential
bias or issue with the sequencing lane / tile of the sample. Unfortunately, 
sequencing lane and tile information is not available for this data

[Github issue showing that SRA does not store this information](https://github.com/ncbi/sra-tools/issues/130)

Particularly this issue was seen in three of the samples
that I have chosen to process at this stage, with the effect
occurring with varying degree. These samples were `SRR8440539`, `SRR8440524`, and `SRR8440538`.


```r
per_base_seq_off <- c('SRR8440539', 'SRR8440524', 'SRR8440538')

df %>% 
  filter(Run %in% per_base_seq_off)
```

We see that two of these sample are from male control samples, 
while one is a female AD sample, so this issue is not contained
to one experimental group.

![](/home/x1/Documents/Weill_Cornell/ANGSD/angsd_project/imgs/per_base_seq_content_fail.png){width=50%}

### **Per Sequence GC Content**

After trimming with `trim_galore`, we can inspect the `multiqc_report.html`
generated from the previous command. Upon inspecting this file, sequence
quality still remains largely high, however we see that one sequence
appears to have a skewed per sequence gc content distribution. We note
that this sample is `SRR8440463`, and will keep an eye on this
sample. 

![](/home/x1/Documents/Weill_Cornell/ANGSD/angsd_project/imgs/per_seq_gc_content.png)

We also see that this same sample has the highest percentage of **duplicate reads**, 
leading us to believe more that this sample might be contaminated by other
sequences that do not model the GC content of the target genome. **We note that
this sample was not of concern earlier when studying the cyclic pattern per base sequence content list**.

![](/home/x1/Documents/Weill_Cornell/ANGSD/angsd_project/imgs/multi_qc_trimmed_gen_stats.png)


```r
df %>% 
  filter(Run == "SRR8440463")
```

We see that this sample is from the experimental condition
(AD), thus this GC content change could be due to disease-specific
modifications in the myeloid cells. 

## **STAR Index Generation**

In order to generate a `STAR` index for later mapping, we will first need to
download the human reference genome, as well 
as the human reference annotation of genes that is compatible with 
this version of the human reference genome.


```bash
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz
wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
```

Now that we have both the reference genome and our features of 
interest, we can generate our `STAR` index using the command below, 
saving our index to the `--genomeDir` directory of `human_genome_38_STAR_index`. 


```bash
STAR --runMode  genomeGenerate \
     --runThreadN 16 \
     --genomeDir human_genome_38_STAR_index \
     --genomeFastaFiles human_genome_38/hg38.fa \
     --sjdbGTFfile human_genome_38/gencode.v36.annotation.gtf \
     --sjdbOverhang 49 
```

**Parameter Choice**: 

* `runThreadN 16` -- Tells `STAR` to use all 16 cores available
  on [my personal machine](https://www.lenovo.com/us/en/desktops-and-all-in-ones/legion-desktops/legion-t-series-towers/Lenovo-Legion-T730-28ICO/p/99LE9700308)
* `sjdbOverhang 49` -- In the `STAR` manual it is recommended to set the `--sjdbOverhang` (splice junction data base overhang) to the **read length - 1**. Since the read lengths associated with my dataset are of length **50bp**, this paramter should be set to **49** for my dataset. 

* Run STAR on a file with the reference



```bash
cell_type="myeloid_cells"
fastq_files=$(ls ${cell_type}/*.fastq.gz)

for fastq_file in $fastq_files; do
  base=$(basename $fastq_file)
  prefix=$(echo $base | egrep -o "SRR[0-9]+")
  
  STAR --runMode  alignReads \
       --runThreadN 12 \
       --genomeDir  human_genome_38_STAR_index \
       --readFilesIn  $fastq_file \
       --readGSM3561843FilesCommand  zcat \
       --outSAMtype  BAM  SortedByCoordinate \
       --outFileNamePrefix  alignments/${prefix}. \
       --outSAMattributes All \
       --outTmpKeep="None"
done 
```


```bash
cell_type="myeloid_cells"
trim_galore_outdir="trimmed_${cell_type}"
trimmed_files=$(ls $trim_galore_outdir/*.fq.gz)
fastq_files=$(ls ${cell_type}/*.fastq.gz)

for fastq_file in $trimmed_files; do
  base=$(basename $fastq_file)
  prefix=$(echo $fastq_file | egrep -o "SRR[0-9]+")
  
  STAR --runMode  alignReads \
       --runThreadN 12 \
       --genomeDir  human_genome_38_STAR_index \
       --readFilesIn  $fastq_file \
       --readFilesCommand  zcat \
       --outSAMtype  BAM  SortedByCoordinate \
       --outFileNamePrefix  trimmed_alignments/${prefix}_trimmed. \
       --outSAMattributes All \
       --outReadsUnmapped="Fastx" \
       --outTmpKeep="None"
done 
```

**Parameter Choice**

Notice the use of `--outSAMattributes All` above, this will add the 
following attributes to our output
  
* **NH**:i:4 -- Number of loci that the read has mapped to.
* **HI**:i:1 -- (Hit Index): Current index of which alignment out
  of NH total alignments we are seeing for the sequencing read
  in question.
* **AS**:i:47	-- Alignment Score
* **nM**:i:1 -- Number of mismatches per paired alignment
* **NM**:i:1 -- Number of mismatches in each mate.
* **MD**:Z:7C42	-- String encoding of the mismatched and deleted reference bases.
* **jM**:B:c,-1	-- Intron motifs for all junctions
* **jI**:B:i,-1 -- Start and End of introns for all junctions (1-based)

`--outSAMtype  BAM  SortedByCoordinate`:  Used so that `samtools sort` is not needed before generating the `.bai` bam-index file with `samtools index` in order to use in most downstream tools. 

`--readFilesCommand  zcat`: Used to maintain compression of `fastq` files
while still being able to run `STAR` over them.

`--outReadsUnmapped="Fastx"`: When aligning reads that may have concerning
quality or length, it might be useful to see which reads (and how many) were
not mapped. This option prints all unmapped reads to a Fasta/FastQ file
with the same prefix for later evaluation.

`--outTmpKeep="None"`: We will not be delving into the inner workings
of STAR, so we will not need to keep any of the temporary intermediate 
files. 

Now we can observe the results of our alignment: 


```bash
samtools view SRR8440443_trimmed.Aligned.sortedByCoord.out.bam | head -n 3
```



```bash
SRR8440443.14857858	16	chr1	10069	1	50M	*	0	0	ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFBFFF<FBB<BB	NH:i:4	HI:i:1	AS:i:47	nM:i:1	NM:i:1	MD:Z:39C10	jM:B:c,-1	jI:B:i,-1
SRR8440443.14631312	16	chr1	10550	0	46M	*	0	0	GTGCAGAGGAGAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBB	NH:i:5	HI:i:1	AS:i:43	nM:i:1	NM:i:1	MD:Z:10C35	jM:B:c,-1	jI:B:i,-1
SRR8440443.16566141	0	chr1	10553	1	49M	*	0	0	CAGAGGAGAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTG	BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NH:i:4	HI:i:1	AS:i:46	nM:i:1	NM:i:1	MD:Z:7C41	jM:B:c,-1	jI:B:i,-1
```


We see that all of the tags mentioned above are included on every alignment In order
to view and later make use our alignment we must index it as well. Note that this
alignment is already sorted as we have included `SortedByCoordinate` in our 
`STAR` command, thus all we must do is run `samtools index` to generate the 
`bai` bam index file.


```bash
bam_files=$(ls *.bam)

for bam in $bam_files ; do
  echo $bam
  samtools index $bam
done
```


### **BAM Quality Control**


```bash
featureCounts \
  -a human_genome_38/gencode.v36.annotation.gtf \
  -o featCounts_genes.txt \
  -T 12 \
  alignments/SRR???????.Aligned.sortedByCoord.out.bam
```

The output below produced from running the above command can also 
be helpful in QCing BAM files, as we can view how many reads 
were successfully assigned alignments. 


```bash
         ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
	  v2.0.1

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 12 BAM files                                     ||
||                           o SRR8440443.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440448.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440463.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440477.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440482.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440501.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440511.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440514.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440518.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440524.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440538.Aligned.sortedByCoord.out.bam       ||
||                           o SRR8440539.Aligned.sortedByCoord.out.bam       ||
||                                                                            ||
||             Output file : featCounts_genes.txt                             ||
||                 Summary : featCounts_genes.txt.summary                     ||
||              Annotation : gencode.v36.annotation.gtf (GTF)                 ||
||      Dir for temp files : ./                                               ||
||                                                                            ||
||                 Threads : 12                                               ||
||                   Level : meta-feature level                               ||
||              Paired-end : no                                               ||
||      Multimapping reads : not counted                                      ||
|| Multi-overlapping reads : not counted                                      ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\============================================================================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file gencode.v36.annotation.gtf ...                        ||
||    Features : 1429877                                                      ||
||    Meta-features : 60660                                                   ||
||    Chromosomes/contigs : 25                                                ||
||                                                                            ||
|| Process BAM file SRR8440443.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 75947707                                             ||
||    Successfully assigned alignments : 15825223 (20.8%)                     ||
||    Running time : 0.08 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440448.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 88271538                                             ||
||    Successfully assigned alignments : 16657112 (18.9%)                     ||
||    Running time : 0.09 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440463.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 79491797                                             ||
||    Successfully assigned alignments : 9882510 (12.4%)                      ||
||    Running time : 0.08 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440477.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 106248338                                            ||
||    Successfully assigned alignments : 16398474 (15.4%)                     ||
||    Running time : 0.11 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440482.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 64754285                                             ||
||    Successfully assigned alignments : 10485796 (16.2%)                     ||
||    Running time : 0.08 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440501.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 92982435                                             ||
||    Successfully assigned alignments : 17984891 (19.3%)                     ||
||    Running time : 0.11 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440511.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 125441957                                            ||
||    Successfully assigned alignments : 14433126 (11.5%)                     ||
||    Running time : 0.14 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440514.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 95716096                                             ||
||    Successfully assigned alignments : 15349065 (16.0%)                     ||
||    Running time : 0.11 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440518.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 84929075                                             ||
||    Successfully assigned alignments : 17075999 (20.1%)                     ||
||    Running time : 0.10 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440524.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 78337110                                             ||
||    Successfully assigned alignments : 9681568 (12.4%)                      ||
||    Running time : 0.09 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440538.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 159236670                                            ||
||    Successfully assigned alignments : 19428974 (12.2%)                     ||
||    Running time : 0.20 minutes                                             ||
||                                                                            ||
|| Process BAM file SRR8440539.Aligned.sortedByCoord.out.bam...               ||
||    Single-end reads are included.                                          ||
||    Total alignments : 52292159                                             ||
||    Successfully assigned alignments : 7792866 (14.9%)                      ||
||    Running time : 0.07 minutes                                             ||
||                                                                            ||
|| Write the final count table.                                               ||
|| Write the read assignment summary.                                         ||
||                                                                            ||
|| Summary of counting results can be found in file "featCounts_genes.txt.su  ||
|| mmary"                                                                     ||
||                                                                            ||
\\============================================================================//
```

We can further inspect why these alignments failed by looking into 
the `featCounts_genes.txt.summary` file. We show the command as 
well as the first entry in the file below. We see that for sample
`SRR8440443`, **15,794,534** mappings were assigned, **43,589,052** alignments
were unassigned due to multi-mapping, **16,357,784** alignments were unassigned
due to "No Features", and **507,488** alignments were unassigned due to ambiguity.


```bash
cat featCounts_genes.txt.summary 
```


```bash
Status	alignments/SRR8440443_trimmed.Aligned.sortedByCoord.out.bam
Assigned	15794534
Unassigned_Unmapped	0
Unassigned_Read_Type	0
Unassigned_Singleton	0
Unassigned_MappingQuality	0
Unassigned_Chimera	0
Unassigned_FragmentLength	0
Unassigned_Duplicate	0
Unassigned_MultiMapping	43589052
Unassigned_Secondary	0
Unassigned_NonSplit	0
Unassigned_NoFeatures	16357784
Unassigned_Overlapping_Length	0
Unassigned_Ambiguity	507488
```

**Samtools Flagstat**

We can also use `samtools flagstat` to determine information 
about our mappings. Below I have written some code to loop
over all of my alignment files, appending the `samtools flagstat` 
output to their respective output files. 


```bash
untrimmed_bam_files=$(ls /home/t730/angsd/alignments | grep -e "SRR[0-9]*.Aligned.*.bam$")
trimmed_bam_files=$(ls /home/t730/angsd/alignments | grep -e "SRR[0-9]*_trimmed.Aligned.*.bam$")

trimmed_out_file="trimmed_bam_quality.txt"
untrimmed_out_file="untrimmed_bam_quality.txt"

rm $trimmed_out_file $untrimmed_out_file

for file in $untrimmed_bam_files ; do   
  echo -e "===================================" >> $untrimmed_out_file
  echo -e "${file}\n" >> $untrimmed_out_file
  samtools flagstat $file >> $untrimmed_out_file
  echo -e "===================================\n" >> $untrimmed_out_file
done

for file in $trimmed_bam_files ; do
  echo -e "===================================" >> $trimmed_out_file
  echo -e "${file}\n" >> $trimmed_out_file
  samtools flagstat $file >> $trimmed_out_file
  echo -e "===================================\n" >> $trimmed_out_file
done
```

This leaves us with entires
in the output files like this


```bash
===================================
SRR8440482_trimmed.Aligned.sortedByCoord.out.bam

64996775 + 0 in total (QC-passed reads + QC-failed reads)
30094645 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
64996775 + 0 mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
===================================
```


**Samtools Coverage**

`samotols coverage` can help us understand the first step of alignment, 
seeing roughly where alighments have been assigned to and plot their 
coverage with the `--histogram` option to print an ascii histogram
representation of the coverage of each chromosome/feature/segment to
the terminal. The **first two** histograms are shown below.


```bash
samtools coverage --histogram SRR8440443_trimmed.Aligned.sortedByCoord.out.bam 
```

![](/home/x1/Documents/Weill_Cornell/ANGSD/angsd_project/imgs/samtools_coverage_example.png)

### **Calculating Unmapped Reads**


```bash
bam_files=$(ls /home/t730/angsd/trimmed_alignments/*.bam)

for bam_file in $bam_files; do
  base=$(basename $bam_file)
  prefix=$(echo $base | egrep -o "SRR[0-9]+")
  
  num_lines=$(cat ${prefix}_trimmed.Unmapped.out.mate1 | wc -l)
  num_total_lines=$(zcat ../trimmed_myeloid_cells/${prefix}_trimmed.fq.gz | wc -l)
  perc_unmapped_reads=$(echo -e "(${num_lines}/${num_total_lines})*100" | bc -l | xargs printf %.2f)
  echo -e "Percentage of unmapped reads for ${prefix}: ${perc_unmapped_reads}%"
done 
```


```bash
Percentage of unmapped reads for SRR8440443: 5.90%
Percentage of unmapped reads for SRR8440448: 4.04%
Percentage of unmapped reads for SRR8440463: 33.79%
Percentage of unmapped reads for SRR8440477: 5.14%
Percentage of unmapped reads for SRR8440482: 6.54%
Percentage of unmapped reads for SRR8440501: 6.98%
Percentage of unmapped reads for SRR8440511: 5.47%
Percentage of unmapped reads for SRR8440514: 4.97%
Percentage of unmapped reads for SRR8440518: 5.02%
Percentage of unmapped reads for SRR8440524: 11.28%
Percentage of unmapped reads for SRR8440538: 7.43%
Percentage of unmapped reads for SRR8440539: 12.34%
```

We see that sample `SRR8440463` showed a very high amount of unmapped reads. 
This is the same sample that showed issues in the **Per Sequence GC Content**
plot before, thus providing us more evident that this sample should be considered
for not being included in the study. 





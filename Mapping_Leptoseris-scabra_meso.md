---
layout: post
title: Mesophotic *Leptoseris cf. scabra* RNAseq data processing
date: '2023-01-05'
categories: Protocols
tags: [RNASeq, Bioinformatics]
---

# Mesophotic *Leptoseris cf. scabra* RNAseq data CBASS heat stress experiment

### Analyzed by Veronica Radice (Barshis Lab), Old Dominion University

## Background
*Field collection*

Fragments of mesophotic colonies of the reef-building coral *Leptoseris cf. scabra* (n=5, 37.3 to 38.7 m) were collected using closed circuit rebreathers at Leone on the island of Tutuila, American Sāmoa in the austral summer (-14.342, -170.789; March 2, 2022). Fragments of each coral were immediately preserved in the field in Zymo DNA/RNA Shield (TField) and the rest of the corals were transported live (~11 minute drive) to the experiment location. 

*Acute heat stress experiment*
An acute heat stress experiment using the CBASS (Coral Bleaching Automated Stress System) was conducted on the same day as the coral collection. Coral fragments were glued to pre-made genotype cards and placed in tanks at ambient temperature (29°C) while the fragments for the initial timepoint (T0) were dark-acclimated for ~30 minutes with an aerator prior to PAM measurements (DIVING PAM II; Walz, Effeltrich, Germany), n=3 measurements per fragment). One fragment of each genotype/species were placed in four tanks with different maximum target temperatures (Control=29°C, Low= 33°C, Medium= 35°C, High=38°C). Temperatures were selected relative to the local maximum monthly mean for American Sāmoa and increased from the control temperature for 3 hours until reaching the target temperature of each tank (profile maintained by Novatech Ice Probes and Bulk Reef Supply 100W Titanium Heaters). Temperature profiles for each tank were controlled via an Arduino with each temperature held (3 hours) followed by a ramp down back to control temperature (~1 hour). Water flow and turnover was maintained in each tank throughout the experiment (SUN JVP Series).

At the end of the maximum temperature hold (T2, after last daylight) 7 hours after the initial timepoint T0, fragments were sampled, photographed, and preserved in Zymo DNA/RNA Shield. The following morning (T3) 18 hours after the initial timepoint T0, PAM measurements were made on each coral fragment prior to being photographed and preserved in DNA/RNA Shield. All preserved coral fragments were frozen at -80°C until subsequent analysis.

--------------------------------------------------------------------------------------------

## Data
Raw data
> /RC/group/rc_barshis_lab/taxonarchive/2022-10-31_BMKData/BMK_DATA_20221021175354_1

Archived files
> /RC/group/rc_barshis_lab/taxonarchive/Leptoseris_scabra/

Data I'm working on
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/2022-10-31_BMKData/

--------------------------------------------------------------------------------------------

## Summary

- RNA was extracted for n=50 *Leptoseris cf. scabra* holobiont samples from 5 individuals (genotypes 1, 4, 6, 9, 10)
- RNA was successfully sequenced for n=46 samples (paired end)

--------------------------------------------------------------------------------------------

## Rename sequencing files

Using /cm/shared/courses/dbarshis/21AdvGenomics/scripts/renamer_advbioinf.py

Already done.

--------------------------------------------------------------------------------------------

# Trimgalore

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/2022-10-31_BMKData/Data/

cp *_Lsca_* /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

```
nano trimgalore_Lsca.sh
```

```
#!/bin/bash -l

#SBATCH -o 2023-01-06_trimgalore_Lsca.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=trimgalore_Lsca

enable_lmod

module load container_env trim_galore

for i in *Lsca*_R1.fq.gz ; do `crun trim_galore --fastqc --paired $i ${i%_R1.fq.gz}_R2.fq.gz`; done
```

```
sbatch trimgalore_Lsca.sh
```


--------------------------------------------------------------------------------------------

# MultiQC

```
enable_lmod
```

```
module load container_env multiqc
```

```
crun multiqc ./
```

Output:
```
  /// MultiQC 🔍 | v1.13

|           multiqc | MultiQC Version v1.14 now available!
|           multiqc | Search path : /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 462/462  
|          cutadapt | Found 92 reports
|            fastqc | Found 92 reports
|           multiqc | Compressing plot data
|           multiqc | Report      : multiqc_report.html
|           multiqc | Data        : multiqc_data
|           multiqc | MultiQC complete
```

--------------------------------------------------------------------------------------------

# Symbiont typing

- based on internal transcribed spacer 2 (ITS2) region, a multicopy genetic marker commonly used to analyse Symbiodinium diversity
- map samples against ITS2 database
- Arif et al. 2014 [https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12869)
- Corrigendum, Arif et al. 2019, with updated ITS2 database file [https://onlinelibrary.wiley.com/doi/10.1111/mec.14956](https://onlinelibrary.wiley.com/doi/10.1111/mec.14956)
- File:  mec14956-sup-0001-files1_corrigendum.fasta

***use Corrigendum Arif ITS2 2019 database***
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/symbiont_typing/

#### Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py mec14956-sup-0001-files1_corrigendum.fasta
```

```
The total number of sequences is 400
The average sequence length is 376
The total number of bases is 150400
The minimum sequence length is 376
The maximum sequence length is 376
The N50 is 376
Median Length = 376
contigs < 150bp = 0
contigs >= 500bp = 0
contigs >= 1000bp = 0
contigs >= 2000bp = 0
```

--------------------------------------------------------------------------------------------

## Map all samples to Arif ITS2 (all clades) 2019 database

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/symbiont_typing/mec14956-sup-0001-files1_corrigendum.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

##### Unzip .gz files
This will extract the compressed file and decompress it.

```
nano gzip_Lsca.sh
```

```
#!/bin/bash -l

#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=gzip_Lsca

gzip -d *_val_*.fq.gz
```

```
sbatch gzip_Lsca.sh
```

#### Build Bowtie index
Build a Bowtie index from a set of DNA sequences (creates 6 files for mapping).
The script bowtie2-build which is part of Bowtie2 requires Python3. 

python --version

python/3.6 (Default)

```
module avail bowtie2
```
bowtie2/2.4 (Default)
samtools/1.10.2

```
nano bowtie2-build.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2-build.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2-build

enable_lmod
module load container_env perl
module load container_env python3
module load bowtie2/2.4

crun.bowtie2 bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/mec14956-sup-0001-files1_corrigendum.fasta Arif_ITS2_corrigendum
```

```
sbatch bowtie2-build.sh
```

### Mapping to Arif_ITS2_corrigendum reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

```
nano mapreads_Arif_ITS2_corrigendum_Lsca.sh
```

Files ending in 
_val_2.fq.gz
_val_1.fq.gz

So use *_val_*.fq

```
#!/bin/bash -l

#SBATCH -o bowtie2_Arif_ITS2_corrigendum_Lsca.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_Arif_ITS2_corrigendum_Lsca

enable_lmod
module load container_env perl
module load container_env python3
module load bowtie2/2.4

for i in *_val_*.fq ; do bowtie2 --rg-id ${i%_*_val_*.fq} \
--rg SM:${i%*_val_*.fq} \
--local -x Arif_ITS2_corrigendum -U $i \
> ${i%*_val_*.fq}_Arif_ITS2_corrigendum.sam -k 5\n; done
```

```
sbatch mapreads_Arif_ITS2_corrigendum_Lsca.sh
```

--------------------------------------------------------------------------------------------

## Count expression - all reads mapped to Arif_ITS2 symbiont database reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

```
nano countexpression_Arif_ITS2_corrigendum.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_Arif_ITS2_corrigendum.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpression_Arif_ITS2_corrigendum_Lsca

enable_lmod
module load container_env python2

crun python2 /cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/*_Arif_ITS2_corrigendum.sam
```

```
sbatch countexpression_Arif_ITS2_corrigendum.sh
```

--------------------------------------------------------------------------------------------

## Merge all _Arif_ITS2_corrigendum_counts.txt files into one big table 

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

- first need to add column to each .txt file with unique sample id
- so we can later identify which sample has which ITS2 sequences

Using regular expressions:

for f in file1 file2 file3; do sed -i "s/$/\t$f/" $f; done

- For each file in the list, this will use sed to append to the end of each line a tab and the filename
- Using the -i flag with sed to perform a replacement in-place, overwriting the file
- Perform a substitution with s/PATTERN/REPLACEMENT/. 
- In this example PATTERN is $, the end of the line, and REPLACEMENT is \t (= a TAB), 
- and $f is the filename, from the loop variable. 
- The s/// command is within double-quotes so that the shell can expand variables
- sed is most practical for pattern substitution and saving in-place. 


```
nano append_filename.sh
```

```
#!/bin/bash -l

#SBATCH -o append_filename.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=append_filename

for f in *_Arif_ITS2_corrigendum_counts.txt; do sed -i "s/$/\t$f/" $f; done
```

```
sbatch append_filename.sh
```

#### Concatenate files

```
cat *_Arif_ITS2_corrigendum_counts.txt > merged_ITS2_counts_Arif-corrigendum_Lsca.txt
```

*Copied file to local folder on computer*
File:  merged_ITS2_counts_Arif-corrigendum_Lsca.xlsx

## Outcome symbiont clade mapping

- 91% reads mapped to I sequences
  - BLAST confirmed 100% match to former Symbiodinium sp. partial rRNA genes and ITS2
- After removing mapping to I sequences
  - 91% mapped to C
- Need to test the few different Cladocopium references
- Then pick the reference with the best mapping rate 

--------------------------------------------------------------------------------------------

## Test symbiont reference transcriptome
### *Cladocopium goreaui* (C1) Davies reference

- Davies et al. (2018). *Symbiodinium functional diversity in the coral Siderastrea siderea is influenced by thermal stress and reef environment, but not ocean acidification.* Frontiers in Marine Science. FMARS-05-00150.
- Assembled and annotated transcriptome for the symbiotic dinoflagellate algae Cladocopium hosted by *Siderastrea siderea* with all host contamination removed. 
- Data were generated using Illumina HiSeq2000 2*100bp reads and assembled using Trinity.
- File:  davies_cladeC_feb.fasta
[https://sites.bu.edu/davieslab/data-code/](https://sites.bu.edu/davieslab/data-code/)
[https://doi.org/10.3389/fmars.2018.00150](https://doi.org/10.3389/fmars.2018.00150)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 65838_davies_cladeC_feb_fixed-clean_suffixed.fasta
```

```
The total number of sequences is 65838
The average sequence length is 1482
The total number of bases is 97581498
The minimum sequence length is 500
The maximum sequence length is 18168
The N50 is 1746
Median Length = 979
contigs < 150bp = 0
contigs >= 500bp = 65838
contigs >= 1000bp = 40840
contigs >= 2000bp = 13246
```

#### Build Bowtie index
Build a Bowtie index from a set of DNA sequences (creates 6 files for mapping).

```
enable_lmod
module load container_env perl
module load container_env python3
module load bowtie2/2.4

crun.bowtie2 bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb_fixed-clean_suffixed.fasta sym_C-goreaui_Davies
```

### Mapping to symbiont *Cladocopium goreaui* Davies transcriptome reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

```
nano mapreads_Cladocopium-goreaui_Davies.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Cladocopium-goreaui_Davies.txt
#SBATCH -n 8
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_Cladocopium-goreaui_Davies

enable_lmod
module load container_env perl
module load container_env python3
module load bowtie2/2.4

for i in *_val_*.fq ; do bowtie2 --rg-id ${i%_*_val_*.fq} \
--rg SM:${i%*_val_*.fq} \
--local -x sym_C-goreaui_Davies -U $i \
> ${i%*_val_*.fq}_sym_C-goreaui_Davies.sam -k 5\n; done
```

```
sbatch mapreads_Cladocopium-goreaui_Davies.sh
```

--------------------------------------------------------------------------------------------

## Count expression - all reads mapped to symbiont *Cladocopium goreaui* Davies transcriptome reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

```
nano countexpression_sym_C-goreaui_Davies.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_sym_C-goreaui_Davies.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpression_sym_C-goreaui_Davies

enable_lmod
module load container_env python2

crun python2 /cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/*_sym_C-goreaui_Davies.sam
```

```
sbatch countexpression_sym_C-goreaui_Davies.sh
```

Outcome:

- avg. 2.75% singly aligned (range 0.6-8.3%)
- avg. 8.83% overall alignment (range 2.1-30.2%)

--------------------------------------------------------------------------------------------

## Test symbiont reference transcriptome 
### *Cladocopium goreaui* (C1) Chen reference 

- Chen et al. (2020). Evidence That Inconsistent Gene Prediction Can Mislead Analysis of Dinoflagellate Genomes. Journal of Phycology.
- predicted protein‐coding genes from published draft Symbiodiniaceae genome of *Cladocopium goreaui* from Liu et al. (2018)
    - Liu et al. (2018). Symbiodinium genomes reveal adaptive evolution of functions related to coral‐dinoflagellate symbiosis. Commun. Biol. 1:95.
        - Symbiodinium goreaui (Clade C, type C1; AIMS-aten-C1-MI-cfu-B2, now AIMS culture collection SCF055-01) is a single-cell monoclonal culture first isolated from the coral Acropora tenuis at Magnetic Island (Queensland, Australia) at 3 m depth; this culture is maintained at the Australian Institute of Marine Science, Townsville, Australia.
        - Illumina HiSeq 2500 platform, 2 × 150 bp reads
    - generated 116.0 Gb (614.6 million reads)
- file:  Cladocopium_goreaui.CDS.fna
[https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947](https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947)
- Chen et al. (2019). Revised genome sequences and annotations of six Symbiodiniaceae taxa.
[https://espace.library.uq.edu.au/view/UQ:8279c9a](https://espace.library.uq.edu.au/view/UQ:8279c9a)

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 39006_Cladocopium_goreaui.CDS.fna
```

```
The total number of sequences is 39006
The average sequence length is 1625
The total number of bases is 63419290
The minimum sequence length is 111
The maximum sequence length is 40317
The N50 is 2388
Median Length = 852
contigs < 150bp = 5
contigs >= 500bp = 31703
contigs >= 1000bp = 21315
contigs >= 2000bp = 9624
```

#### Build Bowtie index
Build a Bowtie index from a set of DNA sequences (creates 6 files for mapping).

```
enable_lmod
module load container_env perl
module load container_env python3
module load bowtie2/2.4

crun.bowtie2 bowtie2-build /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/39006_Cladocopium_goreaui.CDS.fna sym_C-goreaui_Chen
```

### Mapping to symbiont *Cladocopium goreaui* Chen transcriptome reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

```
nano mapreads_Cladocopium-goreaui_Chen.sh
```

```
#!/bin/bash -l

#SBATCH -o bowtie2_Cladocopium-goreaui_Chen.txt
#SBATCH -n 8
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=bowtie2_Cladocopium-goreaui_Chen

enable_lmod
module load container_env perl
module load container_env python3
module load bowtie2/2.4

for i in *_val_*.fq ; do bowtie2 --rg-id ${i%_*_val_*.fq} \
--rg SM:${i%*_val_*.fq} \
--local -x sym_C-goreaui_Chen -U $i \
> ${i%*_val_*.fq}_sym_C-goreaui_Chen.sam -k 5\n; done
```

```
sbatch mapreads_Cladocopium-goreaui_Chen.sh
```

--------------------------------------------------------------------------------------------

## Count expression - all reads mapped to symbiont *Cladocopium goreaui* Chen transcriptome reference

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

```
nano countexpression_sym_C-goreaui_Chen.sh
```

```
#!/bin/bash -l

#SBATCH -o countexpression_sym_C-goreaui_Chen.txt
#SBATCH -n 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=countexpressionsym_C-goreaui_Chen

enable_lmod
module load container_env python2

crun python2 /cm/shared/courses/dbarshis/barshislab/CCourtney/PoloRNASeq/scripts/countxpression_SB_advbioinf.py /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/*_sym_C-goreaui_Chen.sam
```

```
sbatch countexpression_sym_C-goreaui_Chen.sh
```

Outcome:

- avg. 2.84% singly aligned (range 0.1-9.4%)
- avg. 5.85% overall alignment (range 0.5-16.5%)

--------------------------------------------------------------------------------------------

## Decision about symbiont reference 

***Use Cladocopium goreaui Davies transcriptome**

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb_fixed-clean_suffixed.fasta

--------------------------------------------------------------------------------------------

## Coral host de novo genome *Leptoseris scabra* (new 'Wild Genome project' genome)

> /RC/group/rc_barshis_lab/taxonarchive/Leptoseris_scabra/2023-04_WildGenomesLscaAnnotation/

File: PO3051_Leptoseris_scabra.transcript.fasta.gz

## Check assembly details
```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py PO3051_Leptoseris_scabra.transcript.fasta.gz
```

```
The total number of sequences is 225
The average sequence length is 68086
The total number of bases is 15319427
The minimum sequence length is 389
The maximum sequence length is 296486
The N50 is 125671
Median Length = 40318
contigs < 150bp = 0
contigs >= 500bp = 224
contigs >= 1000bp = 223
contigs >= 2000bp = 219
```

cp PO3051_Leptoseris_scabra.transcript.fasta.gz 225_Leptoseris_scabra.transcript.fasta.gz

gzip -d 225_Leptoseris_scabra.transcript.fasta.gz

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py 225_Leptoseris_scabra.transcript.fasta
```

```
The total number of sequences is 35741
The average sequence length is 1586
The total number of bases is 56703834
The minimum sequence length is 59
The maximum sequence length is 57804
The N50 is 2172
Median Length = 990
contigs < 150bp = 509
contigs >= 500bp = 31069
contigs >= 1000bp = 20366
contigs >= 2000bp = 8506
```

**Best practice to rename with contig number**

```
grep -c ">" 225_Leptoseris_scabra.transcript.fasta
```
35741

```
cp 225_Leptoseris_scabra.transcript.fasta 35741_Leptoseris_scabra.transcript.fasta
```

File:  35741_Leptoseris_scabra.transcript.fasta

--------------------------------------------------------------------------------------------

## Concatenate host and symbiont references to make Holobiont reference

Host
> /RC/group/rc_barshis_lab/taxonarchive/Leptoseris_scabra/2023-04_WildGenomesLscaAnnotation/35741_Leptoseris_scabra.transcript.fasta

Symbiont
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb_fixed-clean_suffixed.fasta

```
cat /RC/group/rc_barshis_lab/taxonarchive/Leptoseris_scabra/2023-04_WildGenomesLscaAnnotation/35741_Leptoseris_scabra.transcript.fasta /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/refassembly/65838_davies_cladeC_feb_fixed-clean_suffixed.fasta > Lscabra-transcript_Cgoreaui-Davies_hybridref.fasta
```

## Check assembly details

```
/cm/shared/courses/dbarshis/18AdvBioinf/scripts/avg_cov_len_fasta_advbioinf.py Lscabra-transcript_Cgoreaui-Davies_hybridref.fasta
```

```
The total number of sequences is 101579
The average sequence length is 1518
The total number of bases is 154285332
The minimum sequence length is 59
The maximum sequence length is 57804
The N50 is 1856
Median Length = 2455
contigs < 150bp = 509
contigs >= 500bp = 96907
contigs >= 1000bp = 61206
contigs >= 2000bp = 21752
```

```
mv Lscabra-transcript_Cgoreaui-Davies_hybridref.fasta 101579_Lscabra-transcript_Cgoreaui-Davies_hybridref.fasta
```

--------------------------------------------------------------------------------------------

#### Create a genome index for STAR

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

```
mkdir reference
```

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/reference/

```
nano GenomeGenerate.sh
```

```
#!/bin/bash -l

#SBATCH -o 2023-05-18_Lsca_genomeGenerate.txt
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=Lscagenomegenerate

enable_lmod
module load container_env star

crun STAR --runMode genomeGenerate --runThreadN 4 --genomeDir ./ --genomeFastaFiles /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/101579_Lscabra-transcript_Cgoreaui-Davies_hybridref.fasta --genomeChrBinNbits 16
```

```
sbatch GenomeGenerate.sh
```

*Script failed when only used "STAR ...", instead now need "crun STAR ...". It keeps changing at ODU HPC, very inconsistent & annoying.*


--------------------------------------------------------------------------------------------

# Mapping - STAR Alignment (Spliced Transcripts Alignment to a Reference)
## Map all the samples to concatenated (host plus symbiont) hybrid reference

STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments.
[https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html#:~:text=STAR%20Alignment%20Strategy,Clustering%2C%20stitching%2C%20and%20scoring](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html#:~:text=STAR%20Alignment%20Strategy,Clustering%2C%20stitching%2C%20and%20scoring)

#### Checking version STAR

```
crun STAR --version
```
2.7.10b

*Previous version*
module load star/2.7.3a
module load STAR/2.7.3a

```
module spider STAR
```
--------------------------------------------------------------------------
  star:
--------------------------------------------------------------------------
     Versions:
        star/2.7
        star/2.7.10b

--------------------------------------------------------------------------

```
module spider star/2.7.3a
```
Lmod has detected the following error:  Unable to find:
"star/2.7.3a".


> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/

Note:  I removed --outReadsUnmapped Fastx from script below. 

*If files are gzipped, add*
--readFilesCommand zcat 

```
nano star_Lsca.sh
```

```
#!/bin/bash -l

#SBATCH -o 2023-05-22_LscaSTARMapping.txt
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=STARLsca

enable_lmod
module load container_env star

for i in *_R1_val_1.fq ; do `crun.star STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/reference/ --runThreadN 16 --outSAMattributes All --outSAMattrRGline ID:${i%_R1_val_1.fq} --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn $i ${i%_R1_val_1.fq}_R2_val_2.fq --outFileNamePrefix ${i%_R1_val_1.fq}_Lscabra_CgoreauiDavies`; done
```

```
sbatch star_Lsca.sh
```

***Job ran for 16-17 hours, then "Failed" ExitCode 127, but appears that files are complete? All .bam files ~1-3G***

```
cat 2023-05-22_LglaSTARMapping.txt 
```

```
You are loading a legacy star module, a up to date star will be loaded instead.

Please run star command with "crun.star" prepend to it, for example:

  crun.star STAR

Otherwise you will get "command not found" error.

To avoid this warning, please use following command to load matlab:

  module load container_env star

Lmod has detected the following error: Infinite Load Loop detected for module:
"star/2.7" file: "/usr/share/lmod/lmod/modulefiles/Core/star/2.7.lua"
While processing the following module(s):
    Module fullname  Module Filename
    ---------------  ---------------
    star/2.7         /usr/share/lmod/lmod/modulefiles/Core/star/2.7.lua

/var/spool/slurmd/job10472938/slurm_script: line 13: STAR: command not found
```

```
ls -alh *_CgoreauiDavies*
```
-rwxrwxrwx 1 vradice users  1.1G May 22 12:13 20220302_Meso_Lgla_2_C_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.0K May 22 12:13 20220302_Meso_Lgla_2_C_T2_Lscabra_CgoreauiDaviesLog.final.out
-rwxrwxrwx 1 vradice users  4.3M May 22 12:13 20220302_Meso_Lgla_2_C_T2_Lscabra_CgoreauiDaviesLog.out
-rwxrwxrwx 1 vradice users  2.1K May 22 12:13 20220302_Meso_Lgla_2_C_T2_Lscabra_CgoreauiDaviesLog.progress.out
-rwxrwxrwx 1 vradice users  201K May 22 12:13 20220302_Meso_Lgla_2_C_T2_Lscabra_CgoreauiDaviesSJ.out.tab
-rwxrwxrwx 1 vradice users  2.1G May 22 12:39 20220302_Meso_Lgla_2_C_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.0K May 22 12:39 20220302_Meso_Lgla_2_C_T3_Lscabra_CgoreauiDaviesLog.final.out
-rwxrwxrwx 1 vradice users  4.3M May 22 12:39 20220302_Meso_Lgla_2_C_T3_Lscabra_CgoreauiDaviesLog.out
-rwxrwxrwx 1 vradice users  3.1K May 22 12:39 20220302_Meso_Lgla_2_C_T3_Lscabra_CgoreauiDaviesLog.progress.out
-rwxrwxrwx 1 vradice users  163K May 22 12:39 20220302_Meso_Lgla_2_C_T3_Lscabra_CgoreauiDaviesSJ.out.tab

```
ls -alh *_CgoreauiDaviesAligned.out.bam
```
-rwxrwxrwx 1 vradice users  1.1G May 22 12:13 20220302_Meso_Lgla_2_C_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.1G May 22 12:39 20220302_Meso_Lgla_2_C_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.2G May 22 12:53 20220302_Meso_Lgla_2_H_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.8G May 22 13:24 20220302_Meso_Lgla_2_L_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.5G May 22 13:38 20220302_Meso_Lgla_2_L_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  3.0G May 22 14:04 20220302_Meso_Lgla_2_M_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.9G May 22 14:35 20220302_Meso_Lgla_2_M_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.9G May 22 15:01 20220302_Meso_Lgla_2_T0_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.2G May 22 15:30 20220302_Meso_Lgla_2_TF_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users 1004M May 22 15:47 20220302_Meso_Lgla_3_C_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  988M May 22 16:03 20220302_Meso_Lgla_3_C_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.3G May 22 16:27 20220302_Meso_Lgla_3_H_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  904M May 22 16:46 20220302_Meso_Lgla_3_H_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  928M May 22 17:03 20220302_Meso_Lgla_3_L_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.4G May 22 17:20 20220302_Meso_Lgla_3_L_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.1G May 22 17:36 20220302_Meso_Lgla_3_M_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.3G May 22 17:55 20220302_Meso_Lgla_3_M_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users 1009M May 22 18:22 20220302_Meso_Lgla_3_T0_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.4G May 22 18:43 20220302_Meso_Lgla_3_TF_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  805M May 22 18:58 20220302_Meso_Lgla_5_C_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.4G May 22 19:15 20220302_Meso_Lgla_5_C_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.7G May 22 19:32 20220302_Meso_Lgla_5_H_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  677M May 22 19:48 20220302_Meso_Lgla_5_H_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  774M May 22 20:06 20220302_Meso_Lgla_5_L_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.3G May 22 20:23 20220302_Meso_Lgla_5_L_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.6G May 22 20:51 20220302_Meso_Lgla_5_M_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.4G May 22 21:09 20220302_Meso_Lgla_5_M_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.8G May 22 21:35 20220302_Meso_Lgla_5_T0_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.8G May 22 22:01 20220302_Meso_Lgla_5_TF_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  894M May 22 22:19 20220302_Meso_Lgla_7_C_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.2G May 22 22:39 20220302_Meso_Lgla_7_C_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.5G May 22 23:02 20220302_Meso_Lgla_7_H_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.1G May 22 23:20 20220302_Meso_Lgla_7_H_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.9G May 22 23:50 20220302_Meso_Lgla_7_L_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.3G May 23 00:05 20220302_Meso_Lgla_7_L_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.9G May 23 00:31 20220302_Meso_Lgla_7_M_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.8G May 23 01:04 20220302_Meso_Lgla_7_M_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.7G May 23 01:28 20220302_Meso_Lgla_7_T0_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  2.4G May 23 01:59 20220302_Meso_Lgla_7_TF_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.2G May 23 02:15 20220302_Meso_Lgla_8_C_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  805M May 23 02:31 20220302_Meso_Lgla_8_C_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.5G May 23 02:48 20220302_Meso_Lgla_8_H_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.5G May 23 03:09 20220302_Meso_Lgla_8_L_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.2G May 23 03:26 20220302_Meso_Lgla_8_L_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.5G May 23 03:41 20220302_Meso_Lgla_8_M_T2_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.2G May 23 04:00 20220302_Meso_Lgla_8_M_T3_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  1.9G May 23 04:23 20220302_Meso_Lgla_8_T0_Lscabra_CgoreauiDaviesAligned.out.bam
-rwxrwxrwx 1 vradice users  817M May 23 04:43 20220302_Meso_Lgla_8_TF_Lscabra_CgoreauiDaviesAligned.out.bam


```
enable_lmod
module load container_env multiqc
crun multiqc ./
```

```
  /// MultiQC 🔍 | v1.13

|           multiqc | MultiQC Version v1.14 now available!
|           multiqc | Search path : /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 1300/1300  
|              star | Found 46 reports
|           bowtie2 | Found 3 reports
|          cutadapt | Found 92 reports
|            fastqc | Found 92 reports
|           multiqc | Compressing plot data
|           multiqc | Previous MultiQC output found! Adjusting filenames..
|           multiqc | Use -f or --force to overwrite existing reports instead
|           multiqc | Report      : multiqc_report_1.html
|           multiqc | Data        : multiqc_data_1
|           multiqc | MultiQC complete
```

2.4M May 23 11:05 multiqc_report_1.html

Output:  multiqc_report_Lscabra_Lscabra-Cgoreaui.html


##### Move *Aligned.out.bam files into new directory.

```
mkdir aligned_hybridref
mv *Aligned.out.bam aligned_hybridref
```

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/aligned_hybridref/

--------------------------------------------------------------------------------------------

# Combined STAR deduplication and SALMON quantification
## Count expression - generate counts files

> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/aligned_hybridref/

## SALMON quantification
**What is Salmon?**

Salmon is based on the philosophy of lightweight algorithms, which use the reference transcriptome (in FASTA format) and raw sequencing reads (in FASTQ format) as input, but do not align the full reads. These tools perform both mapping and quantification. Unlike most lightweight and standard alignment/quantification tools, Salmon utilizes sample-specific bias models for transcriptome-wide abundance estimation. Sample-specific bias models are helpful when needing to account for known biases present in RNA-Seq data including:

- GC bias

- positional coverage biases

- sequence biases at 5' and 3' ends of the fragments

- fragment length distribution

- strand-specific methods

If not accounted for, these biases can lead to unacceptable false positive rates in differential expression studies [2]. The Salmon algorithm can learn these sample-specific biases and account for them in the transcript abundance estimates. Salmon is extremely fast at "mapping" reads to the transcriptome and often more accurate than standard approaches [2].

[https://salmon.readthedocs.io/en/latest/salmon.html#](https://salmon.readthedocs.io/en/latest/salmon.html#)

[https://github.com/hbctraining/Intro-to-rnaseq-hpc-salmon-flipped](https://github.com/hbctraining/Intro-to-rnaseq-hpc-salmon-flipped)

[https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/0feca5559bcbde27cbc7634085f07b3624a36f2f/lessons/08_salmon.md](https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/0feca5559bcbde27cbc7634085f07b3624a36f2f/lessons/08_salmon.md)


```
enable_lmod
module load container_env salmon
module avail salmon
```
salmon/1.9.0


CONTINUE HERE - try to dedup & salmon quant multiple files at once


salmon quant to input multiple sequence files.
e.g. salmon quant -l A -i reference -o result -1 {sample1_1.fastq.gz,sample2_1.fastq.gz} -2 {sample1_2.fastq.gz,sample2_2.fastq.gz}

Multiple alignment files

If your alignments for the sample you want to quantify appear in multiple .bam/.sam files, then you can simply provide the Salmon -a parameter with a (space-separated) list of these files. Salmon will automatically read through these one after the other quantifying transcripts using the alignments contained therein. However, it is currently the case that these separate files must (1) all be of the same library type and (2) all be aligned with respect to the same reference (i.e. the @SQ records in the header sections must be identical).

Salmon allows the user to provide a space-separated list of read files to all of it’s options that expect input files (i.e. -r, -1, -2). When the input is paired-end reads, the order of the files in the left and right lists must be the same.

salmon quant -i index -l IU -1 lib_1_1.fq lib_2_1.fq -2 lib_1_2.fq lib_2_2.fq --validateMappings -o out


[https://salmon.readthedocs.io/en/latest/file_formats.html#](https://salmon.readthedocs.io/en/latest/file_formats.html#)


*For loop to run Salmon on multiple files:*

```
#!/bin/bash

for fq in ~/rnaseq/raw_data/*.fq

do

# create a prefix for the output file
samplename=`basename $fq .fq`

# run salmon
salmon quant -i /n/groups/hbctraining/rna-seq_2019_02/reference_data/salmon_index \
 -l A \
 -r $fq \
 -o ${samplename}_salmon \
 --seqBias \
 --useVBOpt \
 --validateMappings

done
```

[https://github.com/hbctraining/Intro-to-rnaseq-hpc-salmon-flipped/blob/main/lessons/09_quasi_alignment_salmon_sbatch.md](https://github.com/hbctraining/Intro-to-rnaseq-hpc-salmon-flipped/blob/main/lessons/09_quasi_alignment_salmon_sbatch.md)

[https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/0feca5559bcbde27cbc7634085f07b3624a36f2f/lessons/08_salmon.md](https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/0feca5559bcbde27cbc7634085f07b3624a36f2f/lessons/08_salmon.md)




--------------------------------------------------------------------------------------------


--------------------------------------------------------------------------------------------

# Archive files to RC/ drive

- bam files

- raw fastq files

- final counts file

CONTINUE HERE

*.fq.gz
*_val_*.fq

copy_rawfastq.sh

#!/bin/bash -l

#SBATCH -o copy_rawfastq.txt
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mail-user=vradice@odu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=copy

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/*.fq.gz /RC/group/rc_barshis_lab/taxonarchive/Leptoseris_scabra/rawfastq/

cp /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Leptoseris_scabra/data_RNAseq/*_val_*.fq /RC/group/rc_barshis_lab/taxonarchive/Leptoseris_scabra/rawfastq/



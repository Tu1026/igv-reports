## Wyatt Lab Specific Notes

### Installation

- git clone this repository on to the server
- cd into the igv-reports directory
- run `pip install -e .` to install the package
- run git pull to update the package

### Preparing the input files

- Do the standard betastasis melting
- Make sure that the melted file has the following columns: Sample_ID CHROM POSITION REF ALT GENE EFFECT Allele_frequency Read_depth (notice only 1 CHROM column)
- Create a new column called "BAM" and fill it with the url to the bam index file for each sample on the server

### Usage

Simply run the command

```
create_report mutations_melted.tsv --flanking 30 --genome hg38 --sequence 2 --begin 3 --end 3 --info-columns Sample_ID CHROM   POSITION        REF     ALT     GENE    EFFECT  Allele_frequency   Read_depth --bam 11 --output <your_output.html>
```

_Notes: The --flanking parameter is the number of bases to show on either side of the mutation. The --genome parameter is the genome build to use. The --sequence parameter is the column number of the sequence (chromosome) column. The --begin parameter is the column number of the start position column. The --end parameter is the column number of the end position column. The --info-columns parameter is a list of the columns to show in the table. The --bam parameter is the column number of the BAM column. The --output parameter is the name of the output file._

#### Once you have the html file

- Press q to blacklist a mutation
- Press w to ignore a mutation
- Press e to whitelist a mutation
- Press b to undo
- Press s to save the whitelist, blacklist, and a curated list

# igv-reports

A Python application to generate self-contained HTML reports that consist of a table of genomic sites or regions and associated IGV views for each site.
The generated HTML page contains all data neccessary for IGV as uuencoded blobs. It can be opened within a web browser as a static page, with no depenency on the original input files.

## Installation

#### Prerequisites

igv-reports **requires Python 3.6** or greater.

As with all Python projects, use of a **virtual environment** is recommended.
Instructions for creating a virtual environment using `conda` follow.

**1.** Install Anaconda from https://docs.anaconda.com/anaconda/

**2.** Create a virtual environment

```bash
conda create -n igvreports python=3.7.1
conda activate igvreports
```

#### Installing igv-reports

```bash
pip install igv-reports
```

igv-reports requires the package _pysam_ version 0.19.1 or greater, which should be installed automatically. However on OSX this sometimes
fails due to missing dependent libraries. This can be fixed following the procedure below, from the pysam
[docs](https://pysam.readthedocs.io/en/latest/installation.html#installation);  
_"The recommended way to install pysam is through conda/bioconda.
This will install pysam from the bioconda channel and automatically makes sure that dependencies are installed.
Also, compilation flags will be set automatically, which will potentially save a lot of trouble on OS X."_

```bash
conda config --add channels r
conda config --add channels bioconda
conda install pysam
```

## Creating a report

A report consists of a table of sites or regions and an associated IGV view for each site. Reports are created with
the command line script `create_report`, or alternatively `python igv_reports/report.py`. Command line arguments are described below.
Although _--tracks_ is optional, a typical report will include at least an alignment track
(BAM or CRAM) file from which the variants were called.

**Arguments:**

- Required
  - **sites** VCF, BED, MAF, BEDPE, or generic tab delimited file of genomic variant sites. Tabix indexed files are supported and strongly recommended for large files.
  - **fasta** Reference fasta file; must be indexed. This argument should be ommited if --genome is used, otherwise it is required.
- The arguments _begin_, _end_, and _sequence_ are required for a generic tab delimited **sites** file.
  - **--begin** INT. Column of start chromosomal position for **sites** file. Used for generic tab delimited input.
  - **--end** INT. Column of end chromosomal position for **sites**. Used for generic tab delimited input.
  - **--sequence** INT. Column of sequence (chromosome) name.
- Optional for generic tab delimited **sites** file

  - **--zero-based** Specify that the position in the **sites** file is 0-based (e.g. UCSC files) rather than 1-based. Default is `false`.

- Optional
  - **--genome** **_New_** An igv.js genome identifier (e.g. hg38). If supplied fasta, ideogram, and the default annotation track for the specified genome will be used.
  - **--tracks** LIST. Space-delimited list of track files, see below for supported formats. If both _tracks_ and _track-config_ are specified _tracks_ will appear first by default.
  - **--track-config** FILE. File containing array of json configuration objects for igv.js tracks. See the [igv.js wiki](https://github.com/igvteam/igv.js/wiki/Tracks-2.0) for more details. This option allows customization of track parameters. When using this option, the track `url` and `indexURL` properties should be set to the paths of the respective files.
  - **--ideogram** FILE. Ideogram file in UCSC cytoIdeo format.
  - **--template** FILE. HTML template file.
  - **--output** FILE. Output file name; default="igvjs_viewer.html".
  - **--info-columns** LIST. Space delimited list of info field names to include in the variant table. If **sites** is a VCF file these are the info ID values. If **sites** is a tab delimited format these are column names.
  - **--info-columns-prefixes** LIST. For VCF based reports only. Space delimited list of prefixes of VCF info field IDs to include in the variant table. Any info field with ID starting with one of the listed values will be included.
  - **--samples** LIST. Space delimited list of sample (i.e. genotypes) names. Used in conjunction with **--sample-columns**.
  - **--sample-columns** LIST. Space delimited list of VCF sample FORMAT field names to include in the variant table. If **--samples** is specified columns will be restricted to those samples, otherwise all samples will be included.
  - **--flanking** INT. Genomic region to include either side of variant; default=1000.
  - **--standalone** Embed all JavaScript referenced via `<script>` tags in the page.
  - **--sort** Applies to alignment tracks only. If specified alignments are initally sorted by the specified option. Supported values include `BASE, STRAND, INSERT_SIZE, MATE_CHR, and NONE`. Default value is `BASE` for single nucleotide variants, `NONE` (no sorting) otherwise. See the igv.js documentation for more information.
  - **--exclude-flags** INT. Value is passed to samtools as "-F" flag. Used to filter alignments. Default value is 1536 which filters alignments marked "duplicate" or "vendor failed". To include all alignments use `--exclude-flags 0`. See [samtools documentation](http://www.htslib.org/doc/samtools-view.html) for more details.
  - **--idlink** URL tempate for information link for VCF ID values. The token $$ will be substituted with the ID value. Example: `--idlink 'https://www.ncbi.nlm.nih.gov/snp/?term=$$'`

**Track file formats:**

Currently supported track file formats are BAM, CRAM, VCF, BED, GFF3, GTF, WIG, and BEDGRAPH. FASTA. BAM, CRAM, and VCF  
files must be indexed. Tabix is supported and it is recommended that all large files be indexed.

## Examples

Data for the examples are available in the github repository [https://github.com/igvteam/igv-reports](https://github.com/igvteam/igv-reports). The repository can be
downloaded as a zip archive here [https://github.com/igvteam/igv-reports/archive/refs/heads/master.zip](https://github.com/igvteam/igv-reports/archive/refs/heads/master.zip).
It is assumed that the examples are run from the root directory of the repository. Output html is written to the [examples directory](https://github.com/igvteam/igv-reports/tree/master/examples)
and can be viewed [here](https://igv.org/igv-reports/examples/1.5.1).

#### NEW (version 1.5.0) - Create a report using a genome identifier: ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_genome.html)\)

```bash
create_report test/data/variants/variants.vcf.gz \
--genome hg38 \
--flanking 1000 \
--info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC \
--tracks test/data/variants/variants.vcf.gz test/data/variants/recalibrated.bam \
--output examples/example_genome.html
```

#### Create a variant report from a VCF file: ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_vcf.html))

```bash

create_report test/data/variants/variants.vcf.gz \
http://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa \
--ideogram test/data/hg38/cytoBandIdeo.txt \
--flanking 1000 \
--info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC \
--samples reads_1_fastq \
--sample-columns DP GQ \
--tracks test/data/variants/variants.vcf.gz test/data/variants/recalibrated.bam test/data/hg38/refGene.txt.gz \
--output examples/example_vcf.html

```

#### Create a variant report with tracks defined in an [igv.js track config json file](https://github.com/igvteam/igv-reports/tree/master/test/data/variants/trackConfigs.json): ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_config.html))

```bash
create_report test/data/variants/variants.vcf.gz \
https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa \
--ideogram test/data/hg38/cytoBandIdeo.txt \
--flanking 1000 \
--info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC \
--track-config test/data/variants/trackConfigs.json \
--output examples/example_config.html
```

#### Create a variant report from a TCGA MAF file: ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_maf.html))

```bash

create_report test/data/variants/tcga_test.maf \
https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta \
--ideogram test/data/hg19/cytoBandIdeo.txt \
--flanking 1000 \
--info-columns Chromosome Start_position End_position Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 dbSNP_RS \
--tracks  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz \
--output examples/example_maf.html

```

#### Create a variant report from a generic tab-delimited file: ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_tab.html))

```bash

create_report test/data/variants/test.maflite.tsv \
https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta \
--sequence 1 --begin 2 --end 3 \
--ideogram test/data/hg19/cytoBandIdeo.txt \
--flanking 1000 \
--info-columns chr start end ref_allele alt_allele \
--tracks https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz \
--output examples/example_tab.html

```

#### NEW (version 1.5.0) - Create a structural variant report from a bedpe file with two locations (BEDPE format): ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_bedpe.html))

```bash

create_report test/data/variants/SKBR3_Sniffles_tra.bedpe \
--genome hg19 \
--flanking 1000 \
--tracks test/data/variants/SKBR3_Sniffles_variants_tra.vcf test/data/variants/SKBR3.ill.bam https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz \
--output examples/example_bedpe.html
```

#### Create a variant report with custom ID link urls: ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_idlink.html))

```bash

create_report test/data/variants/1kg_phase3_sites.vcf.gz \
https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta \
--ideogram test/data/hg19/cytoBandIdeo.txt \
--flanking 1000 \
--tracks test/data/variants/1kg_phase3_sites.vcf.gz test/data/variants/NA12878_lowcoverage.bam https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz \
--idlink 'https://www.ncbi.nlm.nih.gov/snp/?term=$$' \
--output examples/example_idlink.html

```

#### Create a junction report from a splice-junction bed file: ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_junctions.html))

```bash
create_report test/data/junctions/Introns.38.bed \
https://s3.dualstack.us-east-1.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa \
--type junction \
--ideogram test/data/hg38/cytoBandIdeo.txt \
--track-config test/data/junctions/tracks.json \
--info-columns TCGA GTEx variant_name \
--title "Sample A" \
--output examples/example_junctions.html
```

#### Create a report containing wig and bedgraph files

```bash
create_report test/data/wig/regions.bed \
--genome hg19 \
--exclude-flags 512 \
--tracks test/data/wig/ucsc.bedgraph test/data/wig/mixed_step.wig test/data/wig/variable_step.wig \
--output examples/example_wig.html

```

#### Use of `info-columns-prefixes` option. Variant track only, no alignments. ([Example output](https://igv.org/igv-reports/examples/1.5.1/example_ann.html))

```bash
create_report test/data/infofields/consensus.filtered.ann.vcf \
--genome hg19 \
--flanking 1000 \
--info-columns cosmic_gene \
--info-columns-prefixes clinvar \
--tracks test/data/infofields/consensus.filtered.ann.vcf \
--output https://igv.org/igv-reports/examples/1.5.1/example_ann.html
```

#### Use `--exclude-flags` option to include duplicate alignments in report. Default value is 1536 which filters duplicates and vendor-failed reads.

```bash
create_report test/data/dups/dups.bed \
--genome hg19 \
--exclude-flags 512 \
--tracks test/data/dups/dups.bam \
--output examples/example_dups.html

```

#### Converting genomic files to data URIs for use in igv.js

The script `create_datauri` (`python igv_reports/datauri.py`) converts the contents of a file to a data uri for use in igv.js. The datauri will be printed to stdout. _NOTE_ It is not neccessary to run this script explicitly to create a report, it is documented here
for use with stand-alone igv.js.

**Convert a gzipped vcf file to a datauri.**

```bash
create_datauri test/data/variants/variants.vcf.gz

```

**Convert a slice of a local bam file to a datauri.**

```bash
create_datauri --region chr5:474,969-475,009 test/data/variants/recalibrated.bam
```

**Convert a remote bam file to a datauri.**

```bash
create_datauri --region chr5:474,969-475,009 https://1000genomes.s3.amazonaws.com/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
```

## [_Release Notes_](https://github.com/igvteam/igv-reports/releases)

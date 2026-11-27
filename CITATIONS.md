This file documents the software used in this workflow and links them to citation keys in [`CITATIONS.bib`](CITATIONS.bib).
When using it, please ensure that it is up to date and check if all of the listed software was actually used in your version and configuration.
For example, if you use pangenome mapping, you can skip the `bwa` entries, and if you don't use it, you can skip the `vg giraffe` entry.
Also always make sure to name the exact version of the workflow used, giving a tag like `v6.11.4`.
This ensures that readers can trace the exact versions of the packages used by consulting this tag in the git repository or looking up the archived version on Zenodo.

| tool             | citation                            | used for                                                               |
| ---------------- | ----------------------------------- | ---------------------------------------------------------------------- |
| `altair`         | vanderplas_altair_2018              | interactive plots in the report                                        |
| `arriba`         | uhrig_accurate_2021                 | candidate calling of gene fusions                                      |
| `bcftools`       | danecek_twelve_2021                 | read, sort, normalize, write, annotate (binary) variant call format files |
| `bedtools`       | quinlan_bedtools_2010               | handling of various tasks on genomic region files                      |
| `biopython`      | cock_biopython_2009                 | handle fasta files in Python scripts                                   |
| `bwa index`      | li_aligning_2013                    | indexing the reference genome (skip for pangenome mode)                |
| `bwa mem`        | li_aligning_2013                    | mapping DNA reads to the reference genome (skip for pangenome mode)    |
| `CHM-eval`       | li_synthetic-diploid_2018           | synthetic dataset used for benchmarking                                |
| `ColorBrewer`    | brewer_colorbrewer_2013             | good color schemes                                                     |
| `datavzrd`       | wiegand_datavzrd_2025               | interactive tables in the report                                       |
| `delly`          | rausch_delly_2012                   | candidate variant calling structural variants (indels, SVs)            |
| `dna-seq-varlociraptor` | koster_snakemake-workflowsdna-seq-varlociraptor_2026 | this analysis workflow, mention the exact version used |
| `ensembl-vep`    | mclaren_ensembl_2016                | variant annotation with metadata and effect predictions                |
| `fastp`          | chen_fastp_2025                     | sequencing adapter trimming from raw fastq files                       |
| `FastQC`         | noauthor_s-andrewsfastqc_nodate     | quality control of raw fastq files                                     |
| `fgbio`          | homer_fulcrumgenomicsfgbio_2026     | trimming primers from amplicon sequencing                              |
| `freebayes`      | garrison_haplotype-based_2012       | candidate variant calling small variants (SNV, MNVs, indels)           |
| `GATK`           | oconnor_genomics_nodate             | base quality score recalibration in sequence alignment format          |
| `htslib`         | bonfield_htslib_2021                | tabix indexing and bgzip compression                                   |
| `mosdepth`       | pedersen_mosdepth_2018              | coverage calculations to determine covered genome regions              |
| `multiQC`        | ewels_multiqc_2016                  | gathering of quality control output in one report                      |
| `NumPy`          | harris_array_2020                   | numerical calculations in Python scripts on tables and arrays          |
| `Pandas`         | team_pandas-devpandas_2026          | handling of tables in Python scripts                                   |
| `Pandas`         | mckinney_data_2010                  | handling of tables in Python scripts                                   |
| `picardtools`    | noauthor_picard_2019                | marking duplicates in sequence alignment format files                  |
| `pysam`          | noauthor_pysam-developerspysam_nodate | handle (binary) variant call format files in Python scripts          |
| `r-colorramps`   | keitt_colorramps_2024               | color maps for mutational signature analysis                           |
| `r-glmnet`       | friedman_regularization_2010        | generalized linear models for mutational signature analysis            |
| `r-nnls`         | mullen_nnls_2024                    | non-negative least squares for mutational signature analysis           |
| `r-tidyverse`    | wickham_welcome_2019                | good data handling syntax in R scripts                                 |
| `RColorBrewer`   | neuwirth_rcolorbrewer_2022          | good color schemes in R scripts                                        |
| `samtools`       | danecek_twelve_2021                 | read, sort, index, merge (binary) sequence alignment format files      |
| `scikit-learn`   | pedregosa_scikit-learn_2011         | statistical tests in python scripts                                    |
| `siglasso`       | li_using_2020                       | mutational signature analysis                                          |
| `snakemake`      | molder_sustainable_2021             | workflow management system for the whole analysis                      |
| `STAR`           | dobin_star_2013                     | mapping of RNA reads for gene fusion calling                           |
| `statsmodels`    | seabold_statsmodels_2010            | statistical tests in python scripts                                    |
| `umi_tools`      | smith_umi-tools_2017                | parsing of unique molecular identifiers                                |
| `varlociraptor`  | koster_varlociraptor_2020           | probabilistic variant calling and filtering                            |
| `vcf2maf`        | kandoth_mskccvcf2maf_2018           | converting variant call format files into MAF files                    |
| `vega-lite`      | satyanarayan_vega-lite_2017         | interactive plots in the report                                        |
| `vembrane`       | hartmann_insane_2023                | filtering of (binary) variant call format files and converting into tables |
| `vg giraffe`     | siren_pangenomics_2021              | mapping of DNA reads to pangenome                                      |

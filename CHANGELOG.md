# Changelog

### [3.8.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.8.0...v3.8.1) (2022-03-23)


### Bug Fixes

* make the fasta index an explicit input of rule filter_group_regions ([#124](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/124)) ([402a5d9](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/402a5d924855566d51c81d992c7653501f7a634d))

## [3.8.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.7.1...v3.8.0) (2022-03-12)


### Features

* genotyping support, minor fixes, performance improvements and config schema updates, preparation for continuous GIAB benchmarking ([#90](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/90)) ([d39b8ac](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/d39b8ac91e39f8d85e5a6e06d0302c76bf2c2a59))


### Bug Fixes

* ineffective attempt to patch picard syntax usage ([#118](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/118)) ([6d66fd6](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/6d66fd67dbf905cf744e46d001fc6ec00a21f195))
* effective fix of picard syntax in configs and in the mark_duplicates rule ([#119](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/119)) ([13982fc](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/13982fc79ede58b5c3b1e9d31b0b82c974122608))


### [3.7.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.7.0...v3.7.1) (2022-02-23)


### Bug Fixes

* custom-report: Replace indexes of ANN-field by keys ([#113](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/113)) ([c355bc1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/c355bc13a3c3119341842308348a0954fcbc843a))
* fixed tabix indexing for hg19 REVEL scores


## [3.7.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.6.3...v3.7.0) (2022-02-17)


### Features

* Plot REVEL-scores ([#110](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/110)) ([6fb793c](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/6fb793c324babcd3d34ffe5fa11642b982f74ee7))

### Fixes

* refactored primer filtering for having more robust rust script compilation with less external dependencies (#111)

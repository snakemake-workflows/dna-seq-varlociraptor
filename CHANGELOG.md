# Changelog

### [3.13.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.13.0...v3.13.1) (2022-05-12)


### Bug Fixes

* adjust input ([#146](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/146)) ([a70a34b](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/a70a34bd1de37c73bbb3da9bc11e33063c41ab13))

## [3.13.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.12.0...v3.13.0) (2022-05-11)


### Features

* Improve datavzrd tables ([#144](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/144)) ([c6ef41c](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/c6ef41cb861bb085de124c52530cb7d145377625))


### Bug Fixes

* update to latest datavzrd, which fixes various little glitches ([77b01de](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/77b01de2d22f2d7837c293738ecf47d83c044fa3))

## [3.12.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.11.0...v3.12.0) (2022-05-10)


### Features

* render single sample ([#140](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/140)) ([2858ca8](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/2858ca8e0787a19ae3fcc27fb94d757c75c8716e))


### Bug Fixes

* obtain domain for tick plots in datavzrd report from all related columns ([#142](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/142)) ([2e57b12](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/2e57b12cbed2a330801ccda6fe3f45a61e3a6289))

## [3.11.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.10.2...v3.11.0) (2022-05-10)


### Features

* use Varlociraptor 5.0 ([93568d8](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/93568d86d934f07d74e418ecdb027133578bfd23))


### Bug Fixes

* require Varlociraptor 5.0.3 ([22a9b91](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/22a9b91fea59c26a03d3b7dc081a3f3a4b39bfe2))
* use latest bwa wrapper to avoid an issue with mamba pulling an old r-base as picard dependency ([#141](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/141)) ([c533b27](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/c533b2752e9590dc1b5ca6ffea2d1fcb25ce3d93))

### [3.10.2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.10.1...v3.10.2) (2022-05-02)


### Bug Fixes

* defer final output determination to DAG phase ([#137](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/137)) ([b551881](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b55188154b60a35da95cd1fd48f545a33ff7aa9b))

### [3.10.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.10.0...v3.10.1) (2022-04-30)


### Bug Fixes

* Handle multiple samples ([#134](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/134)) ([3a881bf](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/3a881bff2a3e39927ea3e8cd134a32702cb0a59c))

## [3.10.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.9.0...v3.10.0) (2022-04-26)


### Features

* render scenario with yte ([#132](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/132)) ([6ba7339](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/6ba733934fa1ede015cbfdad0680c290377e478d))

## [3.9.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.8.1...v3.9.0) (2022-04-08)


### Features

* use datavzrd for vcf reporting ([#115](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/115)) ([631426f](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/631426f32c236cc4071f19467c2e0ad0db52f5e2))


### Bug Fixes

* enable specifying filter-specific vembrane options beyond the filter expression ([#121](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/121)) ([5833c45](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/5833c45e8f42eea2964c509e98046f29428f6edd))

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

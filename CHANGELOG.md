# Changelog

## [3.18.2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.18.1...v3.18.2) (2022-10-13)


### Bug Fixes

* correctly filter target regions ([#166](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/166)) ([29f64fc](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/29f64fc8ac39a10be27d81d44c0a7bd5eac8ac1f))

## [3.18.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.18.0...v3.18.1) (2022-09-28)


### Bug Fixes

* render group labels in oncoprint gene views as heatmaps ([6041058](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/604105828bdce27b08ae76ef0a1ccb5f38d7f4f1))

## [3.18.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.17.1...v3.18.0) (2022-09-28)


### Features

* update to datavzrd 2.2 ([b1c7963](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b1c79639dfeee529ed7c63065aa9ef4f3a9c8d14))

## [3.17.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.17.0...v3.17.1) (2022-09-27)


### Bug Fixes

* use latest bcftools concat wrapper which fixes a regression that caused it to just pass on the first given bcf. This was leading to missing a **lot** of calls! Make sure to rerun with the following: `-R merge_calls gather_annotated_calls gather_benchmark_calls bcftools_concat gather_calls` if those rules were executed with release 3.16 or 3.17. ([56be4e1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/56be4e12936fa231c37094b2e1daba52dc7780cd))

## [3.17.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.16.0...v3.17.0) (2022-09-23)


### Features

* allow multiple target_regions BED files defined in config.yaml, merging them into one ([#161](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/161)) ([84064ef](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/84064ef5255a6013b84d3702ba42c9ba0abc7bf1))
* more informative file names ([#148](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/148)) ([4b99df9](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/4b99df99e3c3bb3397043a261ee04abc0f00bae6))

## [3.16.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.15.0...v3.16.0) (2022-09-01)


### Features

* add group annotation capabilities and corresponding display thereof in oncoprint views ([#158](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/158)) ([c8aa3f5](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/c8aa3f502f011d25e72aa4e738e67820c85f0dfc))


### Bug Fixes

* update conda envs and wrappers (if necessary) to work properly with strict priorities ([#159](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/159)) ([a5d7be3](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/a5d7be315fff005f1b5f9fb50e7cbd20682d98ca))

## [3.15.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.14.1...v3.15.0) (2022-08-05)


### Features

* enhance oncoprint report ([#156](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/156)) ([c8357fe](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/c8357fed4fb11ffb837652dd7187f8e6e3ee8f8c))

## [3.14.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.14.0...v3.14.1) (2022-07-25)


### Bug Fixes

* datavzrd related issues ([#152](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/152)) ([dfed4e0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/dfed4e03b80a91f9afa62c7f6f870b085762e6fd))

## [3.14.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.13.1...v3.14.0) (2022-05-20)


### Features

* Hide columns in datavzrd report ([#151](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/151)) ([7db7b35](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/7db7b35efb2e12e7ac492066359e75cefad7aebf))
* Remove plotdata from datavzrd template ([#149](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/149)) ([7d184ff](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/7d184ff125bdbe2c442b105332f48d3cc24f8158))


### Bug Fixes

* automatically calculate genomic sites covered for tumor mutational burden calculation ([#130](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/130)) ([7a26dcc](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/7a26dcc2b0324d5cac46d2d8521e67ead0f2211d))

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

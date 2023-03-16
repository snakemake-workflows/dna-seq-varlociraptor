# Changelog

## [3.22.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.21.1...v3.22.0) (2023-03-01)


### Features

* Add datavzrd fields ([#209](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/209)) ([142ead7](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/142ead7b5a987fe48bfdef1239f8e038cb14a36b))


### Bug Fixes

* explicitly use categorical type for chi2 testing ([#210](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/210)) ([079c756](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/079c75627556cc37a0c8c8a4df870170791d033f))
* use varlociraptor 6.0 with its new smart FDR control. This leads to improved filtration results with less false negatives ([#212](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/212)) ([cb00c33](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/cb00c33eb9bb7659c878f9255d7a4915a520686e))

## [3.21.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.21.0...v3.21.1) (2023-02-17)


### Bug Fixes

* fix typo when retrieving candidate filter expression ([dd1e661](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/dd1e66185994b1b79dc0bc3bfa26a5acca636a8a))

## [3.21.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.20.1...v3.21.0) (2023-02-17)


### Features

* Plot simple observations ([#201](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/201)) ([c940ec5](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/c940ec5ef8f9d37d076e2c5422f44a6ce5d171ab))
* Set excluded regions for delly ([#203](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/203)) ([3688ea8](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/3688ea895f660993e4566f681146caee95e7d38f))


### Bug Fixes

* properly handle complex candidate filter expressions including aux files ([#205](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/205)) ([6e353f7](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/6e353f788fd755ec4719d01346a864627bebf734))


### Performance Improvements

* update to varlociraptor 5.10.0 ([#207](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/207)) ([076cf4e](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/076cf4e6806da5c5c9c35f39176edd90be306da0))
* update to varlociraptor 5.9 ([#200](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/200)) ([38135f2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/38135f2a39b09c807b1ba8313a62d9be4d8f5899))

## [3.20.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.20.0...v3.20.1) (2023-02-14)


### Bug Fixes

* properly handle missing candidate filter ([#199](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/199)) ([bdec1a0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/bdec1a0340de290e7ff05ac79162160f6bc39ff8))

## [3.20.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.19.6...v3.20.0) (2023-02-10)


### Features

* add ability to use aux-files for candidate filtering ([80f0ba9](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/80f0ba9677a37963d0c6365d6d366a3c12c658fd))


### Performance Improvements

* update to varlociraptor 5.8.0 ([b15cd85](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b15cd85c49542daf146aeac4afaab3254729a166))

## [3.19.6](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.19.5...v3.19.6) (2023-02-01)


### Bug Fixes

* provide fasta to vep wrapper, update to latest vep wrapper ([814eee3](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/814eee38146d3702666afd1ffb951e00dcf967fa))
* update to latest vep cache wrapper in order to avoid another proxy problem ([8f9d1a2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/8f9d1a2cfa676760db7dbf7271c964a65c2833a2))

## [3.19.5](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.19.4...v3.19.5) (2023-01-30)


### Bug Fixes

* resource download fixes (proxy support for get_vep_cache, curl deployment for download_revel) ([#195](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/195)) ([9e250c6](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/9e250c65eca3f9637fbaf0eda401c80e084ac7db))

## [3.19.4](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.19.3...v3.19.4) (2023-01-23)

### Fixes

* update to Varlociraptor 5.6.2

## [3.19.3](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.19.2...v3.19.3) (2023-01-15)


### Performance Improvements

* update to varlociraptor 5.6.1 ([#189](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/189)) ([4c34c84](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/4c34c84ebe8d3abd590fce23f833f2a617709a39))

## [3.19.2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.19.1...v3.19.2) (2023-01-13)


### Bug Fixes

* reimplement ability to generate testcase; fix environment for tsv_to_excel rule ([#190](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/190)) ([fb2e692](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/fb2e6929af0befb44b4d47fcd8f7b80aa0809d28))

## [3.19.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.19.0...v3.19.1) (2022-12-21)


### Bug Fixes

* typo in datavzrd rule ([19481c2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/19481c23fdaa24f98ee02060deda13bb9c746fa6))
* update to datavzrd 2.11.1 ([a98a80f](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/a98a80f558b40d1f2e6b5307ebcf52741b3f49a8))

## [3.19.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.18.7...v3.19.0) (2022-12-13)


### Features

* Ordering gene oncoprint rows based on annotation label ([#170](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/170)) ([b832e5b](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b832e5bf0a2f37f024f515446ea3771a18c7a771))


### Bug Fixes

* update to varlociraptor 5.4.1 ([#186](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/186)) ([1445902](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/1445902b3d7ab44f2dfd69108e3a844b914c9ce5))

## [3.18.7](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.18.6...v3.18.7) (2022-11-24)


### Bug Fixes

* fail for empty target regions ([#181](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/181)) ([7563aaa](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/7563aaa9d4c5d0f099556525bec2b9633ef32fae))
* Filter offtarget variants ([#180](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/180)) ([3c87042](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/3c87042ff06b493059213f2c8dddf005ea0d3cd8))
* remove clinical significance duplicate ([#182](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/182)) ([33c28c6](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/33c28c6970b0976151181970363153ae5614551d))


### Performance Improvements

* update datavzrd to 2.9 ([#183](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/183)) ([b3e96a2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b3e96a2db09794c02234b3fe5115bf9e314e3672))

## [3.18.6](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.18.5...v3.18.6) (2022-11-15)


### Performance Improvements

* update to varlociraptor 5.3.3 ([#178](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/178)) ([41fe68a](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/41fe68acc975985797cabb8eb9adead10370fc8b))

## [3.18.5](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.18.4...v3.18.5) (2022-11-14)


### Bug Fixes

* fixed output of get_sra rule ([#176](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/176)) ([1cfda3a](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/1cfda3afeb9c9ab6b31667c00fc0c83069f9e648))

## [3.18.4](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.18.3...v3.18.4) (2022-11-11)


### Bug Fixes

* Fix empty regions bed file ([#171](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/171)) ([feac864](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/feac864b0767a30c1fc4faf6e641bb4cd5517900))
* omit software env for some rules when checking between workflow cache ([#168](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/168)) ([d5899da](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/d5899da4163211e90a5d799663f0ad35ae6e7026))

## [3.18.3](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.18.2...v3.18.3) (2022-11-11)


### Bug Fixes

* fix revel score download URL ([#172](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/172)) ([37bc738](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/37bc7381bb93df0ad76d926f13ad889e9a5fc528))

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

# Changelog

## [5.10.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.10.0...v5.10.1) (2024-10-23)


### Bug Fixes

* fix excluding events from the report ([#332](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/332)) ([acdd9e3](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/acdd9e3a87f71e95ca359cf8650ec998f346ff7e))
* more explicity pinning of varlociraptor version ([f66c6a5](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/f66c6a5ef94b5c03878d20cde39ef8fd8688ad6c))
* processing vembrane config ([#330](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/330)) ([1629c6d](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/1629c6d851dcc6336a0e6cb0ac084249c3d427dc))

## [5.10.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.9.0...v5.10.0) (2024-09-26)


### Features

* allow excluding events from the report ([#324](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/324)) ([b5d1966](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b5d196642f7fc56210b8acec4c9379a854e86d88))


### Bug Fixes

* fix testcase rule for remote storage ([#322](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/322)) ([1f2f3b2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/1f2f3b25b41cfabfd42cfa76430ffda22f7247db))
* update to latest datavzrd and vembrane ([#329](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/329)) ([6c5cc97](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/6c5cc97f9c033a7422b41c6108c2a4e608a55b0b))

## [5.9.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.8.4...v5.9.0) (2024-08-26)


### Features

* add maftools output ([#297](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/297)) ([e749260](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/e749260fcfc451a516771d2cbf0c0ac673acb41d))

## [5.8.4](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.8.3...v5.8.4) (2024-08-13)


### Bug Fixes

* make the `aux-files` example use the existing `ANN` field `SYMBOL` ([#320](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/320)) ([0dd0c36](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/0dd0c3673a98231292c1b8abb28c49a039b8dcfb))

## [5.8.3](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.8.2...v5.8.3) (2024-07-24)


### Miscellaneous Chores

* release 5.8.3 ([2fc6a49](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/2fc6a49415ad4650c0873d1036562505a0d6770c))
* update datavzrd wrapper: [#316](https://github.com/snakemake-workflows/dna-seq-varlociraptor/pull/316)

## [5.8.2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.8.1...v5.8.2) (2024-07-18)


### Bug Fixes

* update to latest datavzrd ([b8bd1fc](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b8bd1fcf54ff2e160922b91169c7aa8fcc1aa328))

## [5.8.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.8.0...v5.8.1) (2024-07-17)


### Bug Fixes

* always use the same type for label list in datavzrd config ([#312](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/312)) ([9d20f50](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/9d20f5003da40cfc8414b4877eaab8f2eceaa6dd))
* fix typo in description of variant calls tables ([9d02274](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/9d02274cc722bc55980c7bde303fef42a36ad88e))
* minor fixes ([#311](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/311)) ([a96e071](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/a96e07146aef1a2bc2fd3b85bda3c3ebe48aa0d5))

## [5.8.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.7.0...v5.8.0) (2024-07-01)


### Features

* mutational signatures ([#308](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/308)) ([f037305](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/f0373052750502832c94e6c2bfdc4f1328ef9634))


### Bug Fixes

* fusion sample selection ([#310](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/310)) ([d817454](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/d8174548d1a20dc9bc59150b606caec4d35831c5))

## [5.7.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.6.0...v5.7.0) (2024-06-07)


### Features

* Add improved visualization for clinical significance ([#301](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/301)) ([82eb4dd](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/82eb4dd193db0fd228a286534e91a4eb69476b62))
* add linkouts to non-coding tables ([#302](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/302)) ([15ce609](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/15ce6095fcdbdc70927cb5a17db53f17b221686d))
* support selecting the fdr control mode via the config ([#306](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/306)) ([ac27206](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/ac2720665906b4de0d24e11347d079375e0f60d2))


### Bug Fixes

* download sra ([#304](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/304)) ([f77b6be](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/f77b6be8590b6859c018ed5826e0afa010fba2e9))

## [5.6.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.5.0...v5.6.0) (2024-05-04)


### Features

* reduce fgbio memory usage ([#296](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/296)) ([71611c8](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/71611c80e6712b016f1df5838c4ceb05e059e388))


### Bug Fixes

* fusion candidate filtering ([#299](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/299)) ([6cc658a](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/6cc658aeed60fe09ff7a9ab3bb10b8a94569343a))
* use latest datavzrd wrapper that integrates template rendering ([#300](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/300)) ([0ad923f](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/0ad923fc3f4f5b36da714925194410803786125a))

## [5.5.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.4.1...v5.5.0) (2024-04-09)


### Features

* add ability to store variants in a population database (BCF) ([#288](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/288)) ([e52c29b](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/e52c29bca9ba13795e404c1e65119ff7647576b0))
* add alphamissense linkout ([#294](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/294)) ([a775627](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/a775627fc894c698f587f0bc811d6b6122d6054c))

## [5.4.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.4.0...v5.4.1) (2024-03-28)


### Bug Fixes

* use latest wrapper for obtaining known variants ([#292](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/292)) ([1268dd9](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/1268dd9caf777c4707ea61e47d5267ef509410ae))

## [5.4.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.3.1...v5.4.0) (2024-03-23)


### Features

* support alphamissense and spliceai ([#287](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/287)) ([91a81fa](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/91a81fa0ac6ac31a8ca3a76bc1eb741f2044e734))


### Bug Fixes

* update to latest cutadapt ([#291](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/291)) ([ca73fa7](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/ca73fa7c6e18af5fa3934850ac7128ba9cd9dbcc))
* use latest wrapper for obtaining known variants ([#290](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/290)) ([c8cf4df](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/c8cf4dfbc1f3666db1a58b03afcf190da445fe13))

## [5.3.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.3.0...v5.3.1) (2024-03-12)


### Bug Fixes

* avoid piping if there is no glob pattern in the given fastqs ([#286](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/286)) ([645e5c0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/645e5c06a50d61e3415e9e7aa6b61501264dba5c))
* require Snakemake 8 ([#284](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/284)) ([9a19971](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/9a1997158425bed989e805026213c47250fd8c58))

## [5.3.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.2.0...v5.3.0) (2024-02-27)


### Features

* fusion calling ([#222](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/222)) ([d0b45ba](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/d0b45babb717ee7587fc3dafc359a6bc9494fd94))


### Performance Improvements

* reorganize tests to first only do dry-runs and then run long and short tests in one runner each ([#280](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/280)) ([ee973c6](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/ee973c649ed0e3c6a657e20c2a806247060ab62a))

## [5.2.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.1.0...v5.2.0) (2023-11-30)


### Features

* colored tick plots in datavzrd report ([#279](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/279)) ([757aa0e](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/757aa0e60e275538a3dc8d8fdeb6195fc65f4a40))
* improve and run sra download ci test case ([#276](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/276)) ([074661e](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/074661e3f7f4fa462a1f77c6a3f3d75346a1a6c3))
* map primers by bwa and update datavzrd template ([#275](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/275)) ([12f22e9](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/12f22e9f6c1127546ed9e20add76bdebf8bbec36))


### Bug Fixes

* unique values when extracting samples column for repeated samples (in multiple groups) ([#278](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/278)) ([aaac1b5](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/aaac1b53676879b04f636099975290d1a050fff6))

## [5.1.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.0.2...v5.1.0) (2023-11-10)


### Features

* gene coverage report ([#270](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/270)) ([19ee0eb](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/19ee0eb9fbdb8c1adbb131c649c8c432fcdb8a9e))
* genomenexus link ([#271](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/271)) ([16c4a50](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/16c4a504b2f17ba1cda13142f6d9896bbe4a603f))


### Bug Fixes

* update fastqc wrapper to `v1.10.0` to include memory parsing fix ([#272](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/272)) ([ed8b391](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/ed8b391419e09fee69929eaf7fac63835683dcdf))

## [5.0.2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.0.1...v5.0.2) (2023-10-09)


### Bug Fixes

* remove superfluous import of FTP remote provider ([ae68b90](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/ae68b9079b92ddf6b0ace6888c171381d297ceee))

## [5.0.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v5.0.0...v5.0.1) (2023-08-25)


### Bug Fixes

* general fixes and optimizations ([#263](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/263)) ([cc486d7](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/cc486d704c22779e0dce22da3c5ed17e12b3100c))
* sample has umis ([#265](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/265)) ([914ddfa](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/914ddfa0fa0578896c2d6bcc9863f1273ace04b4))
* switch mem_gb to mem_mb ([#262](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/262)) ([6bda9b5](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/6bda9b5ee49b51bd5a4adfc737c4ee5f3b57018d))

## [5.0.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v4.2.0...v5.0.0) (2023-08-15)


### ⚠ BREAKING CHANGES

* make annotations for candidate variant filtering configurable ([#251](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/251))

### Features

* allow umis in both paired end reads ([#257](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/257)) ([2fdc798](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/2fdc7985127e54c64b24a8ed5dfb11961ccb94d3))
* make annotations for candidate variant filtering configurable ([#251](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/251)) ([b51bd9d](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b51bd9d8c868e67e305ae7de35ec330f264a26da))


### Bug Fixes

* allow using one sample in multiple groups ([#258](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/258)) ([e336f7a](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/e336f7ab64d56c92bf93fe1080e6dd93007e6f11))

## [4.2.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v4.1.0...v4.2.0) (2023-08-10)


### Features

* update envs and wrappers ([#255](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/255)) ([085df20](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/085df20ff1181a9b58d5b34ff281c5fdc2f0b589))

## [4.1.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v4.0.2...v4.1.0) (2023-08-07)


### Features

* Handle trimmed fastqs plus additional fixes ([#252](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/252)) ([13ae724](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/13ae724c58d4d4bfd45e7cc03f7724b9f17d75ef))

## [4.0.2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v4.0.1...v4.0.2) (2023-06-29)


### Bug Fixes

* hide binned vaf ([#250](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/250)) ([fedd817](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/fedd817e081cfe2522a86ed01f338838a9923b4d))
* set id optional ([#247](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/247)) ([b91e3a5](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/b91e3a5e1843a3126cbd5fdc3541829ec0714a89))


### Performance Improvements

* update to varlociraptor 8.3 ([#249](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/249)) ([c81c440](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/c81c4403487c855ad66805c67ae24e083279ff5d))

## [4.0.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v4.0.0...v4.0.1) (2023-06-14)


### Bug Fixes

* update varlociraptor ([#244](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/244)) ([361c9ee](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/361c9ee879a8ba5af06f6de57dd0b548ca667331))


### Performance Improvements

* update to latest datavzrd ([#245](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/245)) ([d4ee98a](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/d4ee98ae724e4ff6cb5144b66612c0a97a2c6ff1))

## [4.0.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.24.3...v4.0.0) (2023-06-12)


### ⚠ BREAKING CHANGES

* use lowercase for all header names in report tables ([#238](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/238))

### Features

* add binned vaf column for sorting by allele frequency ([#240](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/240)) ([a428b30](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/a428b30b2ed268478577f1fe0717927317525873))
* update observation plot to horizontal stacked bars with annotation variable checkbox ([#242](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/242)) ([762d3be](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/762d3be7ef62dcaee86bf01b4a370b371cf8132a))
* use lowercase for all header names in report tables ([#238](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/238)) ([f142d32](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/f142d32f343492788b59832ca85ce83fa3dffa1c))
* use MANE+clinical for selecting representative transcripts ([#235](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/235)) ([0af4ba0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/0af4ba07cfad35afd0d01356776f90e84010f8ba))


### Bug Fixes

* fix mixed up VEP and VEP-cache versions ([#241](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/241)) ([86ccece](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/86ccece456c98c51348682afe6692ae0388b76fb))
* observations plots ([#243](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/243)) ([615d9cc](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/615d9cc24c60fdb0f197181bf7bf759599664b3b))

## [3.24.3](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.24.2...v3.24.3) (2023-05-09)


### Bug Fixes

* use annotated calls ([#234](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/234)) ([0821a01](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/0821a0191efb37198dbd8a07bcd535ae5ae8644b))

## [3.24.2](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.24.1...v3.24.2) (2023-05-05)


### Bug Fixes

* update to latest snpsift, which fixes a bug that before led to loosing some variants ([#232](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/232)) ([5586986](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/558698614e1970630157587411143e2c30b0a4d4))

## [3.24.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.24.0...v3.24.1) (2023-04-28)


### Bug Fixes

* several fixes ([#227](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/227)) ([4c5ab31](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/4c5ab31aa342f30d23411f6bd1a552a2175287f6))
* update datavzrd to 2.18.5 (fixes for empty callsets) ([#229](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/229)) ([f63479e](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/f63479ecb2170a3df10f2d8635ca9d400612b84c))


### Performance Improvements

* use varlociraptor 8.1 ([#228](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/228)) ([e29d9c3](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/e29d9c38d7ee0257512c8b485512ad469d996de6))

## [3.24.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.23.1...v3.24.0) (2023-04-24)


### Features

* Support new varlociraptor 8.0 observations ([#225](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/225)) ([925f736](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/925f73602cd3658d46f76469257ee5aefc99de2d))

## [3.23.1](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.23.0...v3.23.1) (2023-04-17)


### Bug Fixes

* fix vembrane aux ([#223](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/223)) ([245366c](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/245366c0ce55c2bdaff4a4d0cedff4ab56388b4b))

## [3.23.0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/compare/v3.22.0...v3.23.0) (2023-04-13)


### Features

* handle umis ([#213](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/213)) ([61c54e0](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/61c54e0546757e7bdfb0cad8c1aeae31d9731cb6))
* Improved handling of custom annotation fields  ([#221](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/221)) ([540f993](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/540f993fbcc254fcd7f25cf96d27b276c87a077a))


### Bug Fixes

* prenthesize combined filter expressions ([#214](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/214)) ([fae1b0c](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/fae1b0cb79c6e1685a9021efcb8254cf8877e542))
* update mark_duplicates to newest wrapper syntax ([#217](https://github.com/snakemake-workflows/dna-seq-varlociraptor/issues/217)) ([9d3dd30](https://github.com/snakemake-workflows/dna-seq-varlociraptor/commit/9d3dd30d832bb51235a2c11cf81d5a29101463a4))

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

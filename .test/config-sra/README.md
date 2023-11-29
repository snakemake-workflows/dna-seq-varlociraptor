The files used in this test come from the following project:

https://www.ncbi.nlm.nih.gov/bioproject/PRJNA755173

They are basically different (but closely related) strains of yeast that were isolated at the same geolocation. They are described in this manuscript:

> Lee TJ et al., "Extensive sampling of Saccharomyces cerevisiae in Taiwan reveals ecology and evolution of predomesticated lineages.", Genome Res, 2022 May;32(5):864-877, https://doi.org/10.1101/gr.276286.121

A detailed description of samples can be found in Supplementary Table S6 in this document:
https://genome.cshlp.org/content/suppl/2022/05/02/gr.276286.121.DC1/Supplemental_Tables_.xlsx

I basically chose samples with a comparably low coverage via looking through the list at SRA Run Selector, by restricting the `Assay Type` to `wxs` and sorting by the `Bytes` column:
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP334167

Further, I made sure that they:
* had all been collected at the same site (same `Latitude` and `Longitude` values in Supplementary Table S6)
* had all annotations in common, apart from the two variables in Supplementary Table S6 that define the two different groups:
  * `Substrate isolated` defines the group `soil` made up of `PD09A` and `PD12A`
  *`Medium used in isolation` defines the group `medium_L` made up of `PD13B` and `PD12A`
* all belong to the same clonal group, so that they will only differ in very few mutations

The setup reusing `PD12A` in both groups is meant to mimic a tumor-normal setup, where only a panel of normals is available and is reused in every group.
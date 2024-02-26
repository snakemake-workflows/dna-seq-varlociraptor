This workflow generates annotated variant calls that can be viewed in interactive reports, showing all evidence levels provided by Varlociraptor_.
Adapters were removed with Cutadapt_. Reads were mapped with `BWA MEM`_, PCR and optical duplicates were removed with Picard_.
Candidate variant discovery was performed with Freebayes_ and Delly_. Statisticall assessment of variants was conducted with Varlociraptor_.
Fusion resp. variant calling results, sorted by type, event, and impact can be found under Fusion/Variant calls.
The corresponding Varlociraptor_ scenarios, containing the detailed definition of events can be found unter `Variant calling scenarios`_.

.. _Varlociraptor: https://varlociraptor.github.io
.. _BWA MEM: http://bio-bwa.sourceforge.net
.. _Cutadapt: https://cutadapt.readthedocs.io
.. _Picard: https://broadinstitute.github.io/picard
.. _Freebayes: https://github.com/ekg/freebayes
.. _Delly: https://github.com/dellytools/delly

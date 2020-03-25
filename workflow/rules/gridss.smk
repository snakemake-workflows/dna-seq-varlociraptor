jvm_args = f"-Dreference_fasta=resources/genome.fasta -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=true -Dsamjdk.buffer_size=4194304"

rule GridssCollectMetrics:
    input:
        bam="results/recal/{sample}.sorted.bam"
    output:
        insert_size_metrics=temp("tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.insert_size_metrics"),
        histogram=temp("tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.insert_size_histogram.pdf")
    log:
        "log/gridss/collect_metrics/{sample}.log"
    params:
        tmp_dir=temp("tmp/{sample}.sorted.bam.gridss.working"),
        prefix=temp("tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam"),
        maxcoverage=50000,
        metricsrecords=10000000,
        picardoptions=""
    conda:
        "../envs/gridss.yaml"
    shell: """
(gridss gridss.analysis.CollectGridssMetrics \
{jvm_args} \
TMP_DIR={params.tmp_dir} \
ASSUME_SORTED=true \
I={input.bam} \
O={params.prefix} \
THRESHOLD_COVERAGE={params.maxcoverage} \
FILE_EXTENSION=null \
GRIDSS_PROGRAM=null \
PROGRAM=null \
PROGRAM=CollectInsertSizeMetrics \
STOP_AFTER={params.metricsrecords} \
{params.picardoptions}) > {log} 2>&1
        """

rule GridssCollectMetricsAndExtractSVReads:
    input:
        bam="results/recal/{sample}.sorted.bam",
        insert_size_metrics="tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.insert_size_metrics",
    output:
        sv_metrics=temp("tmp/{sample}.sorted.bam.gridss.working/{sample}.sorted.bam.sv_metrics"),
        namedsorted_bam=temp("tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.namedsorted.bam"),
        metrics=temp(multiext("tmp/{sample}.sorted.bam.gridss.working/{sample}.sorted.bam", ".cigar_metrics", ".coverage.blacklist.bed", ".idsv_metrics", ".insert_size_histogram.pdf", ".insert_size_metrics", ".mapq_metrics", ".tag_metrics")),
    log:
        "log/gridss/collect_metrics_and_extract_sv_reads/{sample}.log"
    params:
        dir=temp("tmp/{sample}.sorted.bam.gridss.working"),
        prefix=temp("tmp/{sample}.sorted.bam.gridss.working/{sample}.sorted.bam"),
        tmp_sort=temp("tmp/{sample}.sorted.bam.gridss.working/{sample}.sorted.namedsort"),
        maxcoverage=50000,
        picardoptions="",
    conda:
        "../envs/gridss.yaml"
    threads:
        50
    shell: """
(gridss gridss.CollectGridssMetricsAndExtractSVReads \
{jvm_args} \
TMP_DIR={params.dir} \
ASSUME_SORTED=true \
I={input.bam} \
O={params.prefix} \
THRESHOLD_COVERAGE={params.maxcoverage} \
FILE_EXTENSION=null \
GRIDSS_PROGRAM=null \
GRIDSS_PROGRAM=CollectCigarMetrics \
GRIDSS_PROGRAM=CollectMapqMetrics \
GRIDSS_PROGRAM=CollectTagMetrics \
GRIDSS_PROGRAM=CollectIdsvMetrics \
GRIDSS_PROGRAM=ReportThresholdCoverage \
PROGRAM=null \
PROGRAM=CollectInsertSizeMetrics \
SV_OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0 \
METRICS_OUTPUT={output.sv_metrics} \
INSERT_SIZE_METRICS={input.insert_size_metrics} \
UNMAPPED_READS=false \
MIN_CLIP_LENGTH=5 \
INCLUDE_DUPLICATES=true \
{params.picardoptions} \
| samtools sort \
-n \
-m 100G \
-T {params.tmp_sort} \
-Obam \
-o {output.namedsorted_bam} \
-@ {threads} \
/dev/stdin) > {log} 2>&1
        """


rule GridssComputeSamTags:
    input:
        ref="resources/genome.fasta",
        idx=rules.bwa_index.output,
        namedsorted_bam="tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.namedsorted.bam",
    output:
        coordinate_bam=temp("tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.coordinate.bam")
    log:
        "log/gridss/compute_sam_tags/{sample}.log"
    params:
        working_dir="tmp",
        tmp_sort=temp("tmp/{sample}.sorted.bam.gridss.working/{sample}.sorted.coordinate-tmp"),
        dir=temp("tmp/{sample}.sorted.bam.gridss.working"),
        picardoptions="",
    conda:
        "../envs/gridss.yaml"
    threads:
        3
    shell: """
(gridss gridss.ComputeSamTags \
{jvm_args} \
TMP_DIR={params.dir} \
WORKING_DIR="{params.working_dir}" \
REFERENCE_SEQUENCE={input.ref} \
COMPRESSION_LEVEL=0 \
I={input.namedsorted_bam} \
O=/dev/stdout \
RECALCULATE_SA_SUPPLEMENTARY=true \
SOFTEN_HARD_CLIPS=true \
FIX_MATE_INFORMATION=true \
FIX_DUPLICATE_FLAG=true \
FIX_SA=true \
FIX_MISSING_HARD_CLIP=true \
TAGS=null \
TAGS=NM \
TAGS=SA \
TAGS=R2 \
TAGS=Q2 \
TAGS=MC \
TAGS=MQ \
ASSUME_SORTED=true \
{params.picardoptions} \
| samtools sort \
-m 100G \
-T {params.tmp_sort} \
-Obam \
-o {output.coordinate_bam} \
-@ {threads} \
/dev/stdin) > {log} 2>&1
        """


rule GridssSoftClipsToSplitReads:
    input:
        ref="resources/genome.fasta",
        idx=rules.bwa_index.output,
        coordinate_bam="tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.coordinate.bam"
    output:
        primary_sv=temp("tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.sc2sr.primary.sv.bam"),
        supp_sv=temp("tmp/{sample}.sorted.bam.gridss.working/tmp.{sample}.sorted.bam.sc2sr.supp.sv.bam"),
    log:
        "log/gridss/soft_clips_to_split_reads/{sample}.log"
    params:
        working_dir="tmp",
        picardoptions="",
    conda:
        "../envs/gridss.yaml"
    threads:
        100
    shell: """
(gridss gridss.SoftClipsToSplitReads \
{jvm_args} \
-Xmx20G \
-Dsamjdk.create_index=false \
-Dgridss.gridss.output_to_temp_file=true \
TMP_DIR={params.working_dir} \
WORKING_DIR={params.working_dir} \
REFERENCE_SEQUENCE={input.ref} \
I={input.coordinate_bam} \
O={output.primary_sv} \
OUTPUT_UNORDERED_RECORDS={output.supp_sv} \
WORKER_THREADS={threads} \
{params.picardoptions}) > {log} 2>&1
        """


rule GridssSortSv:
    input:
        supp_sv="{x}.sc2sr.supp.sv.bam"
    output:
        supp_sv=temp("{x}.sc2sr.suppsorted.sv.bam")
    params:
        tmp_sort="{x}.suppsorted.sv-tmp",
    conda:
        "../envs/gridss.yaml"
    threads:
        100
    shell:
        "samtools sort -m 100G -@ {threads} -T {params.tmp_sort} -Obam -o {output.supp_sv} {input.supp_sv}"


rule GridssMergeSupported:
    threads:
        100
    input:
        primary_sv="{p}/tmp.{x}.bam.sc2sr.primary.sv.bam",
        supp_sv="{p}/tmp.{x}.bam.sc2sr.suppsorted.sv.bam",
    output:
        merged=temp("{p}/{x}.bam.sv.bam")
    conda:
        "../envs/gridss.yaml"
    shell:
        "samtools merge -@ {threads} {output.merged} {input.primary_sv} {input.supp_sv}"


rule GridssAssembleBreakends:
    input:
        ref="resources/genome.fasta",
        idx=rules.bwa_index.output,
        bams=lambda wildcards: expand("results/recal/{sample}.sorted.{ending}", sample=get_group_samples(wildcards), ending=["bam", "bam.bai"]),
        svs=lambda wildcards: expand("tmp/{sample}.sorted.bam.gridss.working/{sample}.sorted.bam.sv.{ending}", sample=get_group_samples(wildcards), ending=["bam", "bam.bai"])
    output:
        assembly="tmp/group.{group}.bam.gridss.working/group.{group}.bam"
    log:
        "log/gridss/assemble_breakends/{group}.log"
    params:
        jobindex="0",
        jobnodes="1",
        working_dir="tmp",
        picardoptions="",
        input_args=lambda wildcards: " ".join(expand("INPUT=results/recal/{sample}.sorted.bam", sample=get_group_samples(wildcards)))
    conda:
        "../envs/gridss.yaml"
    threads:
        100
    shell: """
(gridss gridss.AssembleBreakends \
-Dgridss.gridss.output_to_temp_file=true \
{jvm_args} \
-Xmx100g \
JOB_INDEX={params.jobindex} \
JOB_NODES={params.jobnodes} \
TMP_DIR={params.working_dir} \
WORKING_DIR={params.working_dir} \
REFERENCE_SEQUENCE={input.ref} \
WORKER_THREADS={threads} \
O={output.assembly} \
{params.input_args} \
{params.picardoptions}) > {log} 2>&1
        """


rule GridssCollectMetricsGroup:
    input:
        assembly="tmp/group.{group}.bam.gridss.working/group.{group}.bam"
    output:
        idsv_metrics=temp("tmp/group.{group}.bam.gridss.working/group.{group}.bam.idsv_metrics"),
        mapq_metrics=temp("tmp/group.{group}.bam.gridss.working/group.{group}.bam.mapq_metrics"),
        tag_metrics=temp("tmp/group.{group}.bam.gridss.working/group.{group}.bam.tag_metrics"),
    log:
        "log/gridss/collect_metrics_group/{group}.log"
    params:
        prefix="tmp/group.{group}.bam.gridss.working/group.{group}.bam",
        working_dir="tmp/group.{group}.bam.gridss.working",
        maxcoverage=50000,
        picardoptions="",
    conda:
        "../envs/gridss.yaml"
    shell: """
(gridss gridss.analysis.CollectGridssMetrics \
{jvm_args} \
I={input.assembly} \
O={params.prefix} \
THRESHOLD_COVERAGE={params.maxcoverage} \
TMP_DIR={params.working_dir} \
FILE_EXTENSION=null \
GRIDSS_PROGRAM=null \
GRIDSS_PROGRAM=CollectCigarMetrics \
GRIDSS_PROGRAM=CollectMapqMetrics \
GRIDSS_PROGRAM=CollectTagMetrics \
GRIDSS_PROGRAM=CollectIdsvMetrics \
GRIDSS_PROGRAM=ReportThresholdCoverage \
PROGRAM=null \
PROGRAM=CollectAlignmentSummaryMetrics \
{params.picardoptions}) > {log} 2>&1
        """


rule GridssSoftClipsToSplitReadsAssembly:
    input:
        ref="resources/genome.fasta",
        idx=rules.bwa_index.output,
        assembly="tmp/group.{group}.bam.gridss.working/group.{group}.bam",
        idsv_metrics="tmp/group.{group}.bam.gridss.working/group.{group}.bam.idsv_metrics",
        mapq_metrics="tmp/group.{group}.bam.gridss.working/group.{group}.bam.mapq_metrics",
        tag_metrics="tmp/group.{group}.bam.gridss.working/group.{group}.bam.tag_metrics",
    output:
        assembly_primary_sv=temp("tmp/group.{group}.bam.gridss.working/tmp.group.{group}.bam.sc2sr.primary.sv.bam"),
        assembly_supp_sv=temp("tmp/group.{group}.bam.gridss.working/tmp.group.{group}.bam.sc2sr.supp.sv.bam"),
    log:
        "log/gridss/soft_clips_to_split_reads_assembly/{group}.log"
    params:
        working_dir="tmp",
        picardoptions="",
    conda:
        "../envs/gridss.yaml"
    threads:
        100
    shell: """
(gridss gridss.SoftClipsToSplitReads \
{jvm_args} \
-Xmx50G \
-Dgridss.async.buffersize=16 \
-Dsamjdk.create_index=false \
-Dgridss.gridss.output_to_temp_file=true \
TMP_DIR={params.working_dir} \
WORKING_DIR={params.working_dir} \
REFERENCE_SEQUENCE={input.ref} \
WORKER_THREADS={threads} \
I={input.assembly} \
O={output.assembly_primary_sv} \
OUTPUT_UNORDERED_RECORDS={output.assembly_supp_sv} \
REALIGN_ENTIRE_READ=true \
{params.picardoptions}) > {log} 2>&1
        """


rule GridssIdentifyVariants:
    input:
        bams=lambda wildcards: expand("results/recal/{sample}.sorted.bam", sample=get_group_samples(wildcards)),
        assembly_sv="tmp/group.{group}.bam.gridss.working/group.{group}.bam.sv.bam",
        assembly_sv_index="tmp/group.{group}.bam.gridss.working/group.{group}.bam.sv.bam.bai",
        assembly="tmp/group.{group}.bam.gridss.working/group.{group}.bam",
        ref="resources/genome.fasta",
        svs=lambda wildcards: " ".join(expand("tmp/{sample}.sorted.bam.gridss.working/{sample}.sorted.bam.sv.bam", sample=get_group_samples(wildcards))),
        idx=rules.bwa_index.output,
    output:
        unallocated=temp("tmp/group.{group}.vcf.gridss.working/group.{group}.unallocated.vcf")
    log:
        "log/gridss/indentify_variants/{group}.log"
    params:
        input_args=lambda wildcards: " ".join(expand("INPUT=results/recal/{sample}.sorted.bam", sample=get_group_samples(wildcards))),
        working_dir="tmp",
        picardoptions=""
    conda:
        "../envs/gridss.yaml"
    threads:
        100
    shell: """
(gridss gridss.IdentifyVariants \
-Dgridss.output_to_temp_file=true \
{jvm_args} \
-Xmx50g \
TMP_DIR={params.working_dir} \
WORKING_DIR={params.working_dir} \
REFERENCE_SEQUENCE={input.ref} \
WORKER_THREADS={threads} \
{params.input_args} \
ASSEMBLY={input.assembly} \
OUTPUT_VCF={output.unallocated}) > {log} 2>&1
        """

rule GridssAnnotateVariants:
    input:
        bams=lambda wildcards: expand("results/recal/{sample}.sorted.bam", sample=get_group_samples(wildcards)),
        assembly_sv="tmp/group.{group}.bam.gridss.working/group.{group}.bam.sv.bam",
        assembly_sv_index="tmp/group.{group}.bam.gridss.working/group.{group}.bam.sv.bam.bai",
        unallocated="tmp/group.{group}.vcf.gridss.working/group.{group}.unallocated.vcf",
        ref="resources/genome.fasta",
        idx=rules.bwa_index.output,
        svs=lambda wildcards: " ".join(expand("tmp/{sample}.sorted.bam.gridss.working/{sample}.sorted.bam.sv.bam", sample=get_group_samples(wildcards))),
        assembly="tmp/group.{group}.bam.gridss.working/group.{group}.bam",
    output:
        allocated=temp("tmp/group.{group}.vcf.gridss.working/group.{group}.allocated.vcf"),
    log:
        "log/gridss/annotate_variants/{group}.log"
    params:
        input_args=lambda wildcards: " ".join(expand("INPUT=results/recal/{sample}.sorted.bam", sample=get_group_samples(wildcards))),
        working_dir="tmp",
        picardoptions=""
    conda:
        "../envs/gridss.yaml"
    threads:
        7
    shell:        """
(gridss gridss.AnnotateVariants \
-Dgridss.output_to_temp_file=true \
{jvm_args} \
-Xmx50g \
TMP_DIR={params.working_dir} \
WORKING_DIR={params.working_dir} \
REFERENCE_SEQUENCE={input.ref} \
WORKER_THREADS={threads} \
{params.input_args} \
ASSEMBLY={input.assembly} \
INPUT_VCF={input.unallocated} \
OUTPUT_VCF={output.allocated} \
{params.picardoptions}) > {log} 2>&1
        """

rule GridssAnnotateUntemplatedSequence:
    input:
        allocated="tmp/group.{group}.vcf.gridss.working/group.{group}.allocated.vcf",
        ref="resources/genome.fasta",
        idx=rules.bwa_index.output,
    output:
        vcf="results/gridss_vcf/{group}.vcf"
    log:
        "log/gridss/annotate_untemplated_sequences/{group}.log"
    params:
        working_dir="tmp",
        picardoptions=""
    conda:
        "../envs/gridss.yaml"
    threads:
        100
    shell: """
(gridss gridss.AnnotateUntemplatedSequence \
-Dgridss.output_to_temp_file=true \
-Xmx4g \
{jvm_args} \
TMP_DIR={params.working_dir} \
WORKING_DIR={params.working_dir} \
REFERENCE_SEQUENCE={input.ref} \
WORKER_THREADS={threads} \
INPUT={input.allocated} \
OUTPUT={output.vcf} \
{params.picardoptions}) > {log} 2>&1
        """

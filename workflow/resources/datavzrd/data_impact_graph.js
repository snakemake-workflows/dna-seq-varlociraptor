function parse_impact_scores(value) {
    if (!value) return [];
    const records = JSON.parse(value);
    if (!records || records.length === 0) return [];

    // Collect all unique variants across all haplotypes, sorted by genomic position
    const all_variants = [...new Set(
        records.flatMap(r => {
            const trimmed = r.haplotype.replace('c.[', '').replace(']', '');
            return trimmed ? trimmed.split(';') : [];
        })
    )].sort((a, b) =>
        parseInt(a.replace(/[^0-9].*$/, '')) - parseInt(b.replace(/[^0-9].*$/, ''))
    );

    const n_edges = all_variants.length - 1;
    const samples = [...new Set(records.map(r => r.sample))];

    const by_hap_sample = {};
    for (const r of records) {
        by_hap_sample[r.haplotype + '|||' + r.sample] = r;
    }

    const output = [];

    for (const r of records) {
        output.push({ ...r, row_type: 'score_row' });
    }

    const seen_haplotypes = [...new Set(records.map(r => r.haplotype))];

    for (const haplotype of seen_haplotypes) {
        const trimmed = haplotype.replace('c.[', '').replace(']', '');
        const hap_variant_set = new Set(trimmed ? trimmed.split(';') : []);

        const max_frequency = Math.max(...samples.map(s => {
            const r = by_hap_sample[haplotype + '|||' + s];
            return r ? r.frequency : 0;
        }));

        const any_record = by_hap_sample[haplotype + '|||' + samples[0]];

        for (let i = 0; i < n_edges; i++) {
            const variant_from = all_variants[i];
            const variant_to = all_variants[i + 1];
            const from_present = hap_variant_set.has(variant_from);
            const to_present = hap_variant_set.has(variant_to);

            if (!from_present && !to_present) continue;

            let total_reads = 0;
            const per_sample_parts = [];
            for (const sample of samples) {
                const r = by_hap_sample[haplotype + '|||' + sample];
                const count = r ? (parseInt(r.supporting_reads.split(';')[i]) || 0) : 0;
                total_reads += count;
                per_sample_parts.push(`${sample}: ${count}`);
            }

            output.push({
                row_type: 'edge_row',
                haplotype,
                variant_from,
                variant_to,
                from_position: parseInt(variant_from.replace(/[^0-9].*$/, '')),
                to_position: parseInt(variant_to.replace(/[^0-9].*$/, '')),
                total_reads,
                reads_tooltip: per_sample_parts.join(' | '),
                score: any_record ? parseFloat(any_record.score) : 0,
                max_frequency
            });
        }
    }

    return output;
}

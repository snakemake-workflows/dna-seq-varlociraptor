function parse_impact_scores(value) {
    if (!value) {
        return [];
    } else {
        return JSON.parse(value);
    }
}

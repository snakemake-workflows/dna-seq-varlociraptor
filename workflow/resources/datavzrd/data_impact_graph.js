function parse_impact_scores(value) {
    const fixed = value.replace(/'/g, '"').replace(/,(\s*})/g, "$1");
    return JSON.parse(fixed);
}

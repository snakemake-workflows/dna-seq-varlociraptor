function(value) {
    return value
    .split(',')
    .map(part => part.trim())
    .filter(part => part.length > 0)
    .map(score => ({ score }));
}
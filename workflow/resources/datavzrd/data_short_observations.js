function(value) {
  var regex = /([0-9]+)([A|a|R|r]?[E|B|P|S|V|e|b|p|s|v])/g;
  var effects = {
    "E": "Equal",
    "AB": "Barely (Alt)",
    "AP": "Positive (Alt)",
    "AS": "Strong (Alt)",
    "AV": "Very Strong (Alt)",
    "RE": "Equal (Reference)",
    "RB": "Barely (Reference)",
    "RP": "Positive (Reference)",
    "RS": "Strong (Reference)",
    "RV": "Very Strong (Reference)"
}

  var observations = [];
  while ((result = regex.exec(value)) != null) {
    effect = effects[result[2].toUpperCase()]
    var quality = "Low mapping quality";
    if (result[2] == result[2].toUpperCase()) {
        quality = "High mapping quality";
    }
    observations.push({
        "effect": effect,
        "times": parseFloat(result[1]),
        "quality": quality
    })
  }
  return observations
}
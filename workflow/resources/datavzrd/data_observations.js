function(value) {
  var regex = /([0-9]+)([A|a|R|r]?[E|B|P|S|V|e|b|p|s|v])(\.|[0-9]+)(s|p)(\#|\*|\.)(\+|\-|\*)(\>|\<|\*|\!)(\^|\*)(\$|\.)(\*|\.)/g;
  var effects = {
      "E": "Equal",
      "AB": "Barely (Alt)",
      "AP": "Positive (Alt)",
      "AS": "Strong (Alt)",
      "AV": "Very Strong (Alt)",
      "RE": "Equal (Ref)",
      "RB": "Barely (Ref)",
      "RP": "Positive (Ref)",
      "RS": "Strong (Ref)",
      "RV": "Very Strong (Ref)"
  }
  var observations = [];
  var orientation = {
      "*": "unknown",
      ">": "F1R2",
      "<": "F2R1",
      "!": "non standard"
  }

  while ((result = regex.exec(value)) != null) {
    strand = result[6].replace("*", "±")
    edit_distance = result[3].replace(".", "0")
    effect = effects[result[2].toUpperCase()]
    var quality = "low";
    if (result[2] == result[2].toUpperCase()) {
        quality = "high";
    }
    observations.push({
        "strand": strand,
        "effect": effect,
        "times": parseFloat(result[1]),
        "quality": quality,
        "orientation": orientation[result[7]],
        "edit distance": edit_distance
    })
  }
  return observations
}
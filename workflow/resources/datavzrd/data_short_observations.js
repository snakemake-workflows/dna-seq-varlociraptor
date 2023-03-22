function(value) {
  var regex = /([0-9]+)(A|a|R|r)(N|E|B|P|S|V|n|e|b|p|s|v)/g;
  var effects = {
      "N": "None",
      "E": "Equal",
      "B": "Barely",
      "P": "Positive",
      "S": "Strong",
      "V": "Very Strong"
  }
  var allele = {
    "A": "Alt",
    "R": "Ref"
  }

  var observations = [];
  while ((result = regex.exec(value)) != null) {
    effect = effects[result[3].toUpperCase()]
    var quality = "Low mapping quality";
    if (result[3] == result[3].toUpperCase()) {
        quality = "High mapping quality";
    }
    observations.push({
        "allele": allele[result[2].toUpperCase()],
        "effect": effect,
        "times": parseFloat(result[1]),
        "quality": quality
    })
  }
  return observations
}
function(value) {
  var regex = /([0-9]+)(E|B|P|S|V|e|b|p|s|v)/g;
  var effects = {
    "E": "Equal",
    "B": "Barely",
    "P": "Positive",
    "S": "Strong",
    "V": "Very Strong"
  }
  var observations = [];
  let obs = value.split(",")
  for (i=0; i<2; i++) {
    let allele = (i === 0) ? "Reference" : "Alternative";
    while ((result = regex.exec(obs[i])) != null) {
      effect = effects[result[2].toUpperCase()]
      if (effect !== "Equal") {
        effect +=  ` (${allele})`
      }
      var quality = "Low mapping quality";
      if (result[2] === result[2].toUpperCase()) {
          quality = "High mapping quality";
      }
      observations.push({
          "effect": effect,
          "times": parseFloat(result[1]),
          "quality": quality
      })
    }
  }
  return observations
}
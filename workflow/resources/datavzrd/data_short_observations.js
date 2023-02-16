function(value) {
  var regex = /([0-9]+)(N|E|B|P|S|V|n|e|b|p|s|v)/g;
  var effects = {
      "N": "None",
      "E": "Equal",
      "B": "Barely",
      "P": "Positive",
      "S": "Strong",
      "V": "Very Strong"
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
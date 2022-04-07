function(value) {
  var regex = /([0-9]+)(N|E|B|P|S|V|n|e|b|p|s|v)(s|p)(\#|\*|\.)(\+|\-|\*)(\>|\<|\*|\!)(\^|\*)(\$|\.)(\*|\.)/g;
  var effects = {
      "N": "None",
      "E": "Equal",
      "B": "Barely",
      "P": "Positive",
      "S": "Strong",
      "V": "Very Strong"
  }
  var observations = [];
  var orientation = {
      "*": "unknown",
      ">": "F1R2",
      "<": "F2R1",
      "!": "non standard"
  }
  while ((result = regex.exec(value)) != null) {
    strand = result[5].replace("*", "±")
    effect = effects[result[2].toUpperCase()]
    var quality = "Low mapping quality";
    if (result[2] == result[2].toUpperCase()) {
        quality = "High mapping quality";
    }
    observations.push({
        "strand": strand,
        "strand_orientation": strand + ' ' + orientation[result[6]],
        "effect": effect,
        "times": parseFloat(result[1]),
        "quality": quality,
        "orientation": orientation[result[6]]
    })
  }
  return observations
}
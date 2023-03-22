function(value) {
  var regex = /([0-9]+)(A|a|R|r)(N|E|B|P|S|V|n|e|b|p|s|v)(s|p)(\#|\*|\.)(\+|\-|\*)(\>|\<|\*|\!)(\^|\*)(\$|\.)(\*|\.)/g;
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
  var allele = {
    "A": "Alt",
    "R": "Ref"
  }

  while ((result = regex.exec(value)) != null) {
    strand = result[6].replace("*", "±")
    effect = effects[result[3].toUpperCase()]
    var quality = "Low mapping quality";
    if (result[3] == result[3].toUpperCase()) {
        quality = "High mapping quality";
    }
    observations.push({
        "allele": allele[result[2].toUpperCase()],
        "strand": strand,
        "strand_orientation": strand + ' ' + orientation[result[7]],
        "effect": effect,
        "times": parseFloat(result[1]),
        "quality": quality,
        "orientation": orientation[result[7]]
    })
  }
  return observations
}
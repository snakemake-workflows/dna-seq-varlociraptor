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
  var alignment_map = {
    "s": "single end",
    "p": "paired end"
  }
  var locus_map = {
    "#": "alternative",
    "*": "other",
    ".": "no"
  }
  var orientation_map = {
      "*": "unknown",
      ">": "F1R2",
      "<": "F2R1",
      "!": "non standard"
  }
  var indel_ops = {
    "*": "some",
    ".": "no or irrelevant"
  }
  var softclip_map = {
    "$": "yes",
    ".": "no"
  }
  var position_map = {
    "^": "most found",
    "*": "other or irrelevant"
  }
  var observations = [];
  while ((result = regex.exec(value)) != null) {

    var effect = effects[result[2].toUpperCase()]
    var alignment = alignment_map[result[3]]
    var locus = locus_map[result[4]]
    var strand = result[5].replace("*", "±")
    var orientation = orientation_map[result[6]]
    var read_position = position_map[result[7]]
    var softclip = softclip_map[result[8]]
    var indel = indel_ops[result[9]]
    var quality = "Low mapping quality";
    if (result[2] == result[2].toUpperCase()) {
        quality = "High mapping quality";
    }
    observations.push({
        "strand": strand,
        "group": strand + ',' + orientation + ',' + alignment + ',' +  indel + ',' + locus + ',' + softclip + ',' + read_position,
        "info": [
          'Strand: ' + strand,
          'Orientation: ' + orientation,
          'Alignment: ' + alignment,
          'Indel: ' + indel,
          'Locus: ' + locus,
          'Softclip: ' + softclip,
          'Position: ' + read_position
        ],
        "effect": effect,
        "times": parseFloat(result[1]),
        "quality": quality,
        "orientation": orientation,
        "alignment": alignment,
        "locus": locus,
        "indel operations": indel,
        "position": read_position,
        "softclip": softclip,
        "count": 1
    })
  }
  return observations
}
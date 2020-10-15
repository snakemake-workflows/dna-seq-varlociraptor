// customize column_values to display the attributes of your choice to the sidebar
let column_values = ['id', 'position', 'reference', 'alternatives', 'type'];
// customize which parts of the annotation field to display at the sidebar
let ann_values = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

$(document).ready(function () {
    $("html").on('click', '.variant-row', function () {
        let vis_len = $(this).data('vislen');
        if ($(this).data('packed')) {
            for (let t = 1; t <= vis_len; t++) {
                let compressed_specs = $(this).data('vis' + t.toString());
                let unpacker = new jsonm.Unpacker();
                $(this).data('vis' + t.toString(), unpacker.unpack(compressed_specs));
            }
            $(this).data('packed', false);
        }
        let d = $(this).data('description')
        d = d.replace(/, /g,"\",\"");
        d = d.replace("[","[\"");
        d = d.replace("]","\"]")
        let description = JSON.parse(d);

        for (let t = 1; t <= vis_len; t++) {
            let specs = $(this).data('vis' + t.toString());
            console.log(specs);
            specs.data[1].values.forEach(function(a) {
                if (a.row > 0 && Array.isArray(a.flags)) {
                    let flags = {};
                    a.flags.forEach(function(b) {
                        if (b === 1) {
                            flags[b] = "template having multiple segments in sequencing";
                        } else if (b === 2) {
                            flags[b] = "each segment properly aligned according to the aligner";
                        } else if (b === 4) {
                            flags[b] = "segment unmapped";
                        } else if (b === 8) {
                            flags[b] = "next segment in the template unmapped";
                        } else if (b === 16) {
                            flags[b] = "SEQ being reverse complemented";
                        } else if (b === 32) {
                            flags[b] = "SEQ of the next segment in the template being reverse complemented";
                        } else if (b === 64) {
                            flags[b] = "the first segment in the template";
                        } else if (b === 128) {
                            flags[b] = "the last segment in the template";
                        } else if (b === 256) {
                            flags[b] = "secondary alignment";
                        } else if (b === 512) {
                            flags[b] = "not passing filters, such as platform/vendor quality controls";
                        } else if (b === 1024) {
                            flags[b] = "PCR or optical duplicate";
                        } else if (b === 2048) {
                            flags[b] = "vega lite lines";
                        }
                    });
                    a.flags = flags;
                }
            });
            specs.title = 'Sample: ' + $(this).data('vis-sample' + t.toString());
            specs.width = $('#vis' + t.toString()).width() - 40;
            let v = vegaEmbed('#vis' + t.toString(), specs);
        }

        $("#sidebar").empty();
        $.each($(this).data(), function(i, v) {
            if (i !== 'index' && !i.includes("ann") && column_values.includes(i)) {
                $('#sidebar').append('<tr><th class="thead-dark">' + i + '</th><td>' + v + '</td></tr>');
            }
        });
        $("#ann-sidebar").empty();
        let ann_length = $(this).data('annlen');
        let that = this;
        ann_values.forEach(function (x) {
            let name = description[x];
            $('#ann-sidebar').append('<tr>');
            $('#ann-sidebar').append('<th class="thead-dark" style="position: sticky; left:-1px; z-index: 1; background: white">' + name + '</th>');
            for (let j = 1; j <= ann_length; j++) {
                let ix = x + 1;
                let field = 'ann[' + j + '][' + ix + ']';
                let val = $(that).data(field);
                $('#ann-sidebar').append('<td>' + val + '</td>');
            }
            $('#ann-sidebar').append('</tr>');
        });

        $('custom-sidebar').empty()
        var probabilities = []
        $.each($(this).data("info"), function(key, value) {
            if (key.startsWith('PROB_')) {
                var prob_type = key.replace("PROB_", "");
                prob_type = prob_type.charAt(0) + prob_type.substring(1).toLowerCase().replace("_", " ");
                var probability = Math.pow(10, (-1 * parseFloat(value) / 10));
                probabilities.push({
                    "prob_type": prob_type,
                    "prob": probability
                })
            }
        })
        if (probabilities.length > 0) {
            $('#custom-sidebar').append('<div id="Prob">');
            $('#custom-sidebar').append('</div>');
            var probSpec = {
                "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
                "data": {
                    values: probabilities
                },
                "spacing": 20,
                "hconcat": [{
                    "transform": [{
                        "calculate": "datum.prob < 0.01 ? 1 : 0",
                        "as": "alpha"
                    },
                    {
                        "calculate": "datum.prob < 0.01 ? datum.prob : 0.01",
                        "as": "prob"
                    }
                    ],
                    "mark": "circle",
                    "title": null,
                    "encoding": {
                        "y": {
                            "field": "prob_type",
                            "type": "nominal",
                            "title": null,
                            "axis": {
                                "title": false
                            },
                            "sort": "ascending"
                        },
                        "x": {
                            "field": "prob",
                            "type": "quantitative",
                            "title": "Probability",
                            "scale": {
                                "type": "log"
                            }
                        },
                        "color": {
                            "type": "nominal",
                            "field": "alpha",
                            "scale": {
                                "domain": [0, 1],
                                "range": ["#00000000", "#32a852"]
                            },
                            "legend": null
                        },
                        "tooltip": [{
                            "field": "prob_type",
                            "type": "nominal",
                            "title": "Type"
                        }, {
                            "field": "prob",
                            "type": "quantitative",
                            "title": "Probability"
                        }]
                    }
                },
                {
                    "transform": [{
                        "calculate": "datum.prob >= 0.01 ? 1 : 0",
                        "as": "alpha"
                    },
                    {
                        "calculate": "datum.prob >= 0.01 ? datum.prob : 0.01",
                        "as": "prob"
                    }
                    ],
                    "mark": "circle",
                    "encoding": {
                        "y": {
                            "field": "prob_type",
                            "type": "nominal",
                            "axis": null,
                            "sort": "ascending"
                        },
                        "x": {
                            "field": "prob",
                            "type": "quantitative",
                            "title": "Probability",
                            "scale": {
                                "type": "log"
                            }
                        },
                        "color": {
                            "type": "nominal",
                            "field": "alpha",
                            "scale": {
                                "domain": [0, 1],
                                "range": ["#00000000", "#32a852"]
                            },
                            "legend": null
                        },
                        "tooltip": [{
                            "field": "prob_type",
                            "type": "nominal",
                            "title": "Type"
                        }, {
                            "field": "prob",
                            "type": "quantitative",
                            "title": "Probability"
                        }]
                    }
                }
                ],
                "config": {
                    "view": {
                        "stroke": null
                    },
                    "axis": {
                        "grid": false
                    },
                    "concat": {
                        "spacing": 40
                    }
                },
                "encoding": {
                    "tooltip": [{
                        "field": "prob_type",
                        "type": "nominal",
                        "title": "Type"
                    }, {
                        "field": "prob",
                        "type": "quantitative",
                        "title": "Probability"
                    }]
                }
            };
            vegaEmbed('#Prob', probSpec);
        }

        var regex = /([0-9]+)(N|B|P|S|V|n|b|p|s|v)(s|p)(\+|\-|\*)/g;
        var effects = {
            "N": "None",
            "B": "Barely",
            "P": "Positive",
            "S": "Strong",
            "V": "Very Strong"
        }
        var observations = [];
        $.each($(this).data("format")["OBS"], function(sample_name, observation) {
            while ((result = regex.exec(observation)) != null) {
                strand = result[4].replace("*", "±")
                effect = effects[result[2].toUpperCase()]
                var quality = "Low mapping quality";
                if (result[2] == result[2].toUpperCase()) {
                    quality = "High mapping quality";
                }
                observations.push({
                    "sample": sample_name,
                    "strand": strand,
                    "effect": effect,
                    "times": parseFloat(result[1]),
                    "Quality": quality
                })
            }
        })
        if (observations.length > 0) {
            $('#custom-sidebar').append('<div id="Obs">');
            $('#custom-sidebar').append('</div>');
            var obsSpec = {
                "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
                "data": {
                    values: observations
                },
                "mark": "bar",
                "encoding": {
                    "column": {
                        "field": "sample",
                        "type": "ordinal",
                        "title": "Sample"
                    },
                    "row": {
                        "field": "Quality",
                        "type": "ordinal",
                        "title": null
                    },
                    "x": {
                        "field": "strand",
                        "type": "ordinal",
                        "title": "Strand"
                    },
                    "y": {
                        "field": "times",
                        "type": "quantitative",
                        "title": "Counts"
                    },
                    "color": {
                        "field": "effect",
                        "type": "nominal",
                        "title": "Evidence",
                        "sort": ["None", "Barely", "Positive", "Strong", "Very Strong"],
                        "scale": {
                            "domain": ["None", "Barely", "Positive", "Strong", "Very Strong"],
                            "range": ["#999999", "#2ba6cb", "#afdfee", "#ffa3a3", "#ff5555"]
                        }
                    }
                },
                "config": {
                    "axisX": {
                        "labelAngle": 0
                    }
                }
            };
            vegaEmbed('#Obs', obsSpec);
        }

        var binom_dist = []
        $.each($(this).data("format")["AF"]), function(sample_name, allele_frequency) {
            binom_dist = binom_dist.concat(calcBinomDist(parseFloat(allele_frequency), $(this).data("format")["DP"][sample_name]));
        }
        if (binom_dist.length > 0) {
            $('#custom-sidebar').append('<div id="AlleleFreq">');
            $('#custom-sidebar').append('</div>');
            var distSpec = {
                "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
                "data": {
                    values: binom_dist
                },
                "mark": {
                    "type": "line",
                    "interpolate": "linear"
                },
                "encoding": {
                    "column": {
                        "field": "sample",
                        "type": "ordinal",
                        "title": "Sample"
                    },
                    "x": {
                        "field": "frequency",
                        "type": "quantitative",
                        "title": "Frequency"
                    },
                    "y": {
                        "field": "binomProb",
                        "type": "quantitative",
                        "title": "Density"
                    }
                }
            };
            vegaEmbed('#AlleleFreq', distSpec);
        }

        function calcBinomDist(alleleFreq, coverage) {
            var binomValues = [];
            var frequency = 0;
            var k = Math.round(alleleFreq * coverage);
            while (frequency <= 1) {
                var value = BinomialCoefficient(coverage, k) * Math.pow(frequency, k) * Math.pow((1 - frequency), (coverage - k))
                binomValues.push({
                    "frequency": frequency,
                    "binomProb": value,
                });
                frequency += 0.01;
            }
            return binomValues;
        }

        function BinomialCoefficient(n, k)  {
            if (k + k > n) { k = n - k }
            if (k < 0) { return 0 }
            else {
              var result = 1
              for (i=0;i<k;) {
                result = (result * (n-i)) / ++i
              }
              return result
            }
          }
    })
})
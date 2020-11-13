// customize column_values to display the attributes of your choice to the sidebar
let column_values = ['id', 'position', 'reference', 'alternatives', 'type'];
// customize which parts of the annotation field to display at the sidebar
let ann_values = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

let score_thresholds = {}
score_thresholds["PolyPhen"] = { "Benign": [0, 0.149], "Possibly Damaging": [0.15, 0.849], "Probably Damaging": [0.85, 1] }
score_thresholds["SIFT"] = {"Benign": [0, 0.05], "Damaging": [0.051, 1] }

let score_scales = {}
score_scales["SIFT"] = { "colors": ["#2ba6cb", "#ff5555"], "entries": ["Benign", "Damaging"] }
score_scales["PolyPhen"] = { "colors": ["#2ba6cb", "#ffa3a3", "#ff5555"], "entries": ["Benign", "Possibly Damaging", "Probably Damaging"] }

$(document).ready(function () {



    $("html").on('click', '.variant-row', function () {
        let vis_len = $(this).data('vislen');
        if ($(this).data('packed')) {
            for (let t = 1; t <= vis_len; t++) {
                let compressed_specs = $(this).data('vis' + t.toString());
                let unpacker = new jsonm.Unpacker();
                unpacker.setMaxDictSize(100000);
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


        $('#ann-sidebar').append('<tr>');
        $('#ann-sidebar').append('<th class="thead-dark" style="position: sticky; left:-1px; z-index: 1; background: white">Linkouts</th>');
        var sift_scores = []
        var polyphen_scores = []
        for (let j = 1; j <= ann_length; j++) {
            var transcript = $(that).data('ann[' + j + '][7]')
            sift_score = $(that).data('ann[' + j + '][35]')
            if (sift_score != "") {
                sift_score = sift_score.split("(")[1].slice(0, -1)
                sift_scores = parse_score(sift_scores, sift_score, "SIFT", transcript)
            }

            polyphen_score = $(that).data('ann[' + j + '][36]')
            if (polyphen_score != "") {
                polyphen_score = polyphen_score.split("(")[1].slice(0, -1)
                polyphen_scores = parse_score(polyphen_scores, polyphen_score, "PolyPhen", transcript)
            }

            gene_field = 'ann[' + j + '][4]'
            ensembl_field = 'ann[' + j + '][5]'
            gene = $(that).data(gene_field)
            ensembl_id = $(that).data(ensembl_field)
            $('#ann-sidebar').append('<td id="Linkout' + j +'">')
            $('#Linkout'+ j).append('<div id="Div' + j +'" class="dropdown show">')
            $('#Div'+ j).append('<a class="btn btn-secondary dropdown-toggle" href="#" role="button" id="dropdownMenuLink" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">Select source</a>')
            $('#Div'+ j).append('<div id="Button' + j + '" div class="dropdown-menu" aria-labelledby="dropdownMenuLink">')
            $('#Button'+ j).append('<a class="dropdown-item" href="https://clinicaltrials.gov/ct2/results?cond=&term=' + gene + '&cntry=&state=&city=&dist=" target="_blank">ClinicalTrials</a><br>');
            $('#Button'+ j).append('<a class="dropdown-item" href="https://www.cbioportal.org/ln?q=' + gene + '" target="_blank">cBioPortal</a><br>');
            $('#Button'+ j).append('<a class="dropdown-item" href="https://oncokb.org/gene/' + gene + '" target="_blank">OncoKB</a><br>');
            $('#Button'+ j).append('<a class="dropdown-item" href="https://www.ensembl.org/homo_sapiens/Gene/Summary?db=core;g='+ ensembl_id +'" target="_blank">Ensembl</a>');
            $('#Button'+ j).append('<a class="dropdown-item" href="https://whi.color.com/gene/'+ ensembl_id +'" target="_blank">FLOSSIES</a>');
            $('#Button'+ j).append('<a class="dropdown-item" href="https://varsome.com/gene/'+ gene +'" target="_blank">VarSome</a>');
        }
        $('#ann-sidebar').append('</tr>');


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
                    "mark": {
                        "type": "circle",
                        "size": 45,
                    },
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
                                "range": ["#00000000", "#4682b4"]
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

        var regex = /([0-9]+)(N|B|P|S|V|n|b|p|s|v)(s|p)(\+|\-|\*)(\>|\<|\*|\!)/g;
        var effects = {
            "N": "None",
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
                    "strand_orientation": strand + ' ' + orientation[result[5]], 
                    "effect": effect,
                    "times": parseFloat(result[1]),
                    "quality": quality,
                    "orientation": orientation[result[5]] 
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
                        "field": "quality",
                        "type": "ordinal",
                        "title": null
                    },
                    "x": {
                        "field": "strand_orientation",
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
                    },
                    "tooltip": [ {
                        "field": "strand",
                        "type": "nominative",
                        "title": "Strand"
                    }, {
                        "field": "orientation",
                        "type": "nominal",
                        "title": "Orientation"
                    }, {
                        "field": "times",
                        "type": "quantitative",
                        "title": "Counts"
                    }]
                },
                "config": {
                    "axisX": {
                        "labelAngle": -45
                    }
                }
            };
            vegaEmbed('#Obs', obsSpec);
        }

        var binom_dist = []
        af = $(this).data("format")["DP"]
        $.each($(this).data("format")["AF"], function(sample_name, allele_frequency) {
            binom_dist = binom_dist.concat(calcBinomDist(parseFloat(allele_frequency), af[sample_name][0], sample_name));
        })
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
                    },
                    "tooltip": [{
                        "field": "binomProb",
                        "type": "quantitative",
                        "title": "Density"
                    }, {
                        "field": "frequency",
                        "type": "quantitative",
                        "title": "Frequency"
                    }]
                }
            };
            vegaEmbed('#AlleleFreq', distSpec);
        }

        if (sift_scores.length > 0) {
            $('#custom-sidebar').append('<div id="SiftScores">');
            $('#custom-sidebar').append('</div>');
            plotScores("SIFT", sift_scores, "SiftScores")
        }

        if (polyphen_scores.length > 0) {
            $('#custom-sidebar').append('<div id="PolyPhenScores">');
            $('#custom-sidebar').append('</div>');
            plotScores("PolyPhen", polyphen_scores, "PolyPhenScores")
        }

        function calcBinomDist(alleleFreq, coverage, sample_name) {
            var binomValues = [];
            var frequency = 0;
            var k = Math.round(alleleFreq * coverage);
            while (frequency <= 1) {
                var value = math.combinations(coverage, k) * Math.pow(frequency, k) * Math.pow((1 - frequency), (coverage - k))
                binomValues.push({
                    "frequency": frequency,
                    "binomProb": value,
                    "sample": sample_name
                });
                frequency += 0.01;
            }
            return binomValues;
        }


        function parse_score(scores, variant_value, score_type, transcript) {
            if (score_type == "SIFT") {
                score = 1 - parseFloat(variant_value)
            } else {
                score = parseFloat(variant_value)
            }
            var effect_bins = build_effect_bins(score_type, score)
            for (bin of effect_bins) {
                var effect_idx = score_scales[score_type]["entries"].indexOf(bin.effect)
                scores.push({
                    "name": score_type,
                    "transcript": transcript,
                    "size": bin.size,
                    "effect": bin.effect,
                    "effect_idx": effect_idx
                })
            }
            return scores
        }



        function build_effect_bins(name, score) {
            var effects = score_thresholds[name]
            var effect_bins = []
            if ((Object.keys(effects).length == 1) || score == 0) {
                var effect = extract_effect(name, score)
                effect_bins.push({ "effect": effect, "size": score })
            } else {
                var score_greater_zero = (score, left_boundary, right_boundary) => {
                    if (score < right_boundary) {
                        return (left_boundary > 0) ? (score - left_boundary) : score
                    } else {
                        return (left_boundary > 0) ? (right_boundary - left_boundary) : right_boundary
                    }
                }
                var score_less_zero = (score, left_boundary, right_boundary) => {
                    if (score > left_boundary) {
                        return (right_boundary < 0) ? (score - right_boundary) : score
                    } else {
                        return (right_boundary < 0) ? (left_boundary - right_boundary) : left_boundary
                    }
                }
                for (var effect in effects) {
                    var effect_boundarys = effects[effect]
                    if (score > 0) {
                        if ((effect_boundarys[1] > 0) && (score > effect_boundarys[0])) {
                            var size = score_greater_zero(score, effect_boundarys[0], effect_boundarys[1])
                            effect_bins.push({ "effect": effect, "size": size })
                        }
                    } else if (score < 0) {
                        if ((effect_boundarys[0] < 0) && (score < effect_boundarys[1])) {
                            var size = score_less_zero(score, effect_boundarys[0], effect_boundarys[1])
                            effect_bins.push({ "effect": effect, "size": size })
                        }
                    } else {
                        alert(`Value error: ${score}`)
                    }
                }
            }
            return effect_bins
        }

        function extract_effect(name, value) {
            if (!(name in score_thresholds)) {
                return "Effect not found"
            }
            var effects = score_thresholds[name]
            for (var key in effects) {
                if ((value >= effects[key][0]) && (value <= effects[key][1])) {
                    return key
                }
            }
            return "Score out of range"
        }

        function plotScores(score_type, scores, cell_id) {
            var x_title = {"SIFT": "1-score", "PolyPhen": "Score"}
            var ScoreSpec = {
                "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
                "data": {
                    values: scores
                },
                "mark": "bar",
                "encoding": {
                    "y": {
                        "field": "transcript",
                        "type": "nominal",
                        "title": "Transcript"
                    },
                    "x": {
                        "field": "size",
                        "type": "quantitative",
                        "scale": {"domain": [0, 1]},
                        "title": x_title[score_type],
                    },
                    "color": {
                        "field": "effect",
                        "type": "nominal",
                        "title": "Effect",
                        "sort": score_scales[score_type]["entries"],
                        "scale": {
                            "domain": score_scales[score_type]["entries"],
                            "range": score_scales[score_type]["colors"]
                        }
                    },
                    "order": { "field": "effect_idx", "type": "quantitative" },
                    "tooltip": null
                },
                "config": {
                    "axisX": {
                        "labelAngle": 90
                    }
                }
            };
            vegaEmbed(`#${cell_id}`, ScoreSpec)
        }
    })
})

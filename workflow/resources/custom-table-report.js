let score_thresholds = {}
score_thresholds["PolyPhen"] = { "Benign": [0, 0.149], "Possibly Damaging": [0.15, 0.849], "Probably Damaging": [0.85, 1] }
score_thresholds["SIFT"] = {"Benign": [0, 0.05], "Damaging": [0.051, 1] }
score_thresholds["REVEL"] = {"Likely Benign": [0, 0.5], "Likely Malign": [0.501, 1]}

let score_scales = {}
score_scales["SIFT"] = { "colors": ["#2ba6cb", "#ff5555"], "entries": ["Benign", "Damaging"] }
score_scales["PolyPhen"] = { "colors": ["#2ba6cb", "#ffa3a3", "#ff5555"], "entries": ["Benign", "Possibly Damaging", "Probably Damaging"] }
score_scales["REVEL"] = { "colors": ["#2ba6cb", "#ff5555"], "entries": ["Likely Benign", "Likely Malign"] }

$(document).ready(function () {
    $("html").on('click', '.variant-row', function () {
        let ann_length = $(this).data('annlen');
        let that = this;
        $('#ann-sidebar').append('<tr>');
        $('#ann-sidebar').append('<th class="thead-dark" style="position: sticky; left:-1px; z-index: 1; background: white">Linkouts</th>');
        var sift_scores = []
        var polyphen_scores = []
        var revel_scores = []
        for (let j = 1; j <= ann_length; j++) {
            var transcript_index = ANN_DESCRIPTION.indexOf("Feature")+1
            var transcript = $(that).data('ann[' + j + '][' + transcript_index + ']')
            sift_index = ANN_DESCRIPTION.indexOf("SIFT")+1
            sift_score = $(that).data('ann[' + j + ']['+ sift_index + ']')
            if (sift_score != "" && sift_score !== undefined) {
                sift_score = sift_score.split("(")[1].slice(0, -1)
                sift_scores = parse_score(sift_scores, sift_score, "SIFT", transcript)
            }

            var polyphen_index = ANN_DESCRIPTION.indexOf("PolyPhen")+1
            var polyphen_score = $(that).data('ann[' + j + ']['+ polyphen_index + ']')
            if (polyphen_score != "" && polyphen_score !== undefined) {
                polyphen_score = polyphen_score.split("(")[1].slice(0, -1)
                polyphen_scores = parse_score(polyphen_scores, polyphen_score, "PolyPhen", transcript)
            }

            var revel_index = ANN_DESCRIPTION.indexOf("REVEL")+1
            if (revel_index != 0) {
                var revel_score = $(that).data('ann[' + j + ']['+ revel_index + ']')
                if (revel_score != "" && revel_score !== undefined) {
                    revel_scores = parse_score(revel_scores, revel_score, "REVEL", transcript)
                }
            }

            var linkout_button = true
            $('#ann-sidebar').append('<td id="Linkout' + j +'">')
            $('#Linkout'+ j).append('<div id="Div' + j +'" class="dropdown show">')
            var gene_index = ANN_DESCRIPTION.indexOf("SYMBOL")+1
            gene_field = 'ann[' + j + '][' + gene_index + ']'
            gene = $(that).data(gene_field)
            if (gene !== undefined) {
                linkout_button = create_linkout_button(linkout_button, j);
                $('#Button'+ j).append('<a class="dropdown-item" href="https://clinicaltrials.gov/ct2/results?cond=&term=' + gene + '&cntry=&state=&city=&dist=" target="_blank">ClinicalTrials</a>');
                $('#Button'+ j).append('<a class="dropdown-item" href="https://www.cbioportal.org/ln?q=' + gene + '" target="_blank">cBioPortal</a>');
                $('#Button'+ j).append('<a class="dropdown-item" href="https://oncokb.org/gene/' + gene + '" target="_blank">OncoKB</a>');
                $('#Button'+ j).append('<a class="dropdown-item" href="https://varsome.com/gene/'+ gene +'" target="_blank">VarSome</a>');
                $('#Button'+ j).append('<a class="dropdown-item" href="https://pandrugs.org/#!/query/?genes='+ gene +'" target="_blank">PanDrugs</a>');
                $('#Button'+ j).append('<a class="dropdown-item" href="https://www.mycancergenome.org/content/gene/'+ gene +'" target="_blank">My Cancer Genome</a>');
                $('#Button'+ j).append('<a class="dropdown-item" href="https://search.cancervariants.org/#'+ gene +'" target="_blank">MetaKB</a>');
                $('#Button'+ j).append('<a class="dropdown-item" href="https://www.proteinatlas.org/search/'+ gene +'" target="_blank">The Human Protein Atlas</a>');
                $('#Button'+ j).append('<a class="dropdown-item" href="https://varseak.bio/search.php?gene='+ gene +'" target="_blank">varSEAK Variant Table</a>');
                var hgvsc_index = ANN_DESCRIPTION.indexOf("HGVSc")+1
                transcript_hgvsc_field = 'ann[' + j + '][' + hgvsc_index + ']'
                transcript_hgvsc = $(that).data(transcript_hgvsc_field)
                if (transcript_hgvsc !== undefined) {
                    transcript = transcript_hgvsc.split(":")[0]
                    hgvs = transcript_hgvsc.split(":")[1]
                    $('#Button'+ j).append('<a class="dropdown-item" href="https://varseak.bio/ssp.php?gene=' + gene + '&transcript='+ transcript +'&variant=&hgvs='+ hgvs +'" target="_blank">varSEAK Splice Site Prediction</a>');
                }
            }
            var ensembl_idx = ANN_DESCRIPTION.indexOf("Gene")+1
            ensembl_field = 'ann[' + j + '][' + ensembl_idx + ']'
            ensembl_id = $(that).data(ensembl_field)
            if (ensembl_id !== undefined) {
                create_linkout_button(linkout_button, j)
            }
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
        $.each($(this).data("format")["OBS"], function(sample_name, observation) {
            while ((result = regex.exec(observation)) != null) {
                strand = result[5].replace("*", "±")
                effect = effects[result[2].toUpperCase()]
                var quality = "Low mapping quality";
                if (result[2] == result[2].toUpperCase()) {
                    quality = "High mapping quality";
                }
                observations.push({
                    "sample": sample_name,
                    "strand": strand,
                    "strand_orientation": strand + ' ' + orientation[result[6]],
                    "effect": effect,
                    "times": parseFloat(result[1]),
                    "quality": quality,
                    "orientation": orientation[result[6]]
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
        $.each($(this).data("format")["AFD"], function(sample_name, afd_string) {
            if (afd_string !== undefined) {
                var allele_frequency_distribution = afd_string.split(',')
                for (let dist_idx in allele_frequency_distribution) {
                    var values = allele_frequency_distribution[dist_idx].split('=')
                    var frequency = values[0]
                    var posterior_density = Math.pow(10, (-values[1]/10))
                    binom_dist.push({
                        "frequency": frequency,
                        "binomProb": posterior_density,
                        "sample": sample_name
                    });
                }
            }
        })
        if (binom_dist.length > 0) {
            $('#custom-sidebar').append('<div id="AlleleFreq">');
            $('#custom-sidebar').append('</div>');
            var distSpec = {
                "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
                "data": {
                    values: binom_dist
                },
                "mark": "point",
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
                        "title": "Density",
                        "scale": {"type": "log"}
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

        if (revel_scores.length > 0) {
            $('#custom-sidebar').append('<div id="REVELScores">');
            $('#custom-sidebar').append('</div>');
            plotScores("REVEL", revel_scores, "REVELScores")
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
            var x_title = {"SIFT": "1-score", "PolyPhen": "Score", "REVEL": "Score"}
            var ScoreSpec = {
                "$schema": "https://vega.github.io/schema/vega-lite/v3.json",
                "title": score_type,
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

        function create_linkout_button(linkout_button, j) {
            if (linkout_button == true) {
                $('#Div'+ j).append('<a class="btn btn-secondary dropdown-toggle" href="#" role="button" id="dropdownMenuLink" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">Select source</a>')
                $('#Div'+ j).append('<div id="Button' + j + '" div class="dropdown-menu" aria-labelledby="dropdownMenuLink">')
                return false;
            }
            return true;
        }
    })
})

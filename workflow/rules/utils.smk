import os
from pathlib import Path
import random
import numpy as np
from typing import List, Dict
from collections import defaultdict


###############
### Classes ###
###############

# Implementation of switch/case in Python
# https://code.activestate.com/recipes/410692/

class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    
    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False


#########################
### General Functions ###
#########################

def flatten(
        list_of_lists: List
) -> List:
        """
        Flatten a list of lists recursively

        https://stackoverflow.com/a/53778278

        :param list_of_lists: A list of lists
        :return result: A string that has been flattened from a list of lists
        """

        result = list()
        for i in list_of_lists:
                if isinstance(i, list):
                        result.extend(flatten(i))
                else:
                        result.append(str(i))
        return result


def partial_format(
    string: str,
    **wildcards: str
) -> str:
    """
    Format string while ignoring missing string argument

    https://stackoverflow.com/a/43526674

    :param string: String to be formatted
    :param **wildcards: Optional strings to be formatted
    :return string: Partially formatted string
    """

    class SafeDict(dict):
        def __missing__(self, key):
            return '{' + key + '}'
    replacements = SafeDict(**wildcards)
    return str(string).format_map(replacements)


def partial_expand(
    string: str,
    **wildcards: str
) -> str:
    """Expand the wildcards in `string` from the ones present in wildcards

    Based on function `partially_expand` from:
    https://github.com/snakemake/snakemake/blob/07f73006dcee1e6ce6c271355fff83705f188de4/snakemake/rules.py#L160-L173

    This is done by replacing all wildcard delimiters by `{{` or `}}`
    that are not in `wildcards.keys()`.
    """
    s = str(string).replace("{", "{{").replace("}", "}}")
    for key in wildcards.keys():
          s = s.replace("{{{{{}}}}}".format(key), "{{{}}}".format(key))
    return expand(s, **wildcards)



def to_dict(
    keys: List,
    values: List
) -> Dict:
    """Convert two lists (keys and values) to a dictionary

    Even though it is a simple function, it makes the code more readable.
    """
    return dict(zip(keys, values))


###########################
### Snakemake Functions ###
###########################

def expand_ext(path, ext):
    if isinstance(path, list):
        l = [partial_expand(Path(p).with_suffix('.{ext}'), ext=ext) for p in path]
        l = list(map(list, zip(*l)))
        return flatten(l)
    elif isinstance(path, str):
        return partial_expand(Path(path).with_suffix('.{ext}'), ext=ext)


def ext_dict(path):
    # Exract extensions
    ext = [''.join(Path(p).suffixes) for p in path]
    # Fix extensions to be used as keys
    keys = [e.replace(".gz", "").split(".")[-1] for e in ext]
    return to_dict(keys, path)


def get_tmp(
    name: str = None,
    local: bool = True
) -> str:
    """
    Returns path to temporary folder

    :param name: name of folder (randomly generated if not provided)
    :param local: should the temp folder be local (/tmp/) or not?
    :return path: a string with path to temp folder
    """

    if not name:
        name = "temp_{0:06d}".format(random.randrange(1e6))

    if False and '__default__' in cluster_config:
        path = Path('/scratch/$PBS_JOBID')
        path = Path('/dev/shm')
    else:
        if local:
            path = Path('/tmp')
        else:
            path = Path('temp/large_files')
            path.mkdir(parents=True, exist_ok=True)

    return(str(path / name))


#############
### Rules ###
#############

rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi"
    log:
        "logs/bcf-index/{prefix}.log"
    wrapper:
        "0.67.0/bio/bcftools/index"


rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bai"
    log:
        "logs/bam-index/{prefix}.log"
    wrapper:
        "0.67.0/bio/samtools/index"


rule cram_index:
    input:
        "{prefix}.cram"
    output:
        "{prefix}.crai"
    log:
        "logs/cram-index/{prefix}.log"
    wrapper:
        "0.67.0/bio/samtools/index"


rule tabix_known_variants:
    input:
        "resources/{prefix}.{format}.gz"
    output:
        "resources/{prefix}.{format}.gz.tbi"
    log:
        "logs/tabix/{prefix}.{format}.log"
    params:
        get_tabix_params
    cache: True
    wrapper:
        "0.59.2/bio/tabix"


rule testcase:
    input:
        obs=get_group_observations,
        scenario="results/scenarios/{group}.yaml"
    output:
        directory("resources/testcases/{group}.{caller}/{locus}")
    log:
        "logs/varlociraptor/testcase/{group}.{caller}.{locus}.log"
    params:
        obs=lambda w, input: ["{}={}".format(s, f) for s, f in zip(get_group_aliases(w), input.obs)],
        parent=lambda w, output: os.path.dirname(output[0])
    threads: workflow.cores
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants --testcase-prefix {output} --testcase-locus {wildcards.locus} "
        "generic --obs {params.obs} "
        "--scenario {input.scenario} 2> {log}"

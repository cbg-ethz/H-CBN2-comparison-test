import pandas as pd

# OPTIONS
configfile: "config.yaml"
DATASETS = pd.DataFrame(config["datasets"], index=['1', '2'], columns=['file'])

def get_input(wildcards):
    return DATASETS.loc[wildcards.file, 'file']

rule all:
    input:
        "empirical_distribution.rds"

localrules: split

# 1. Merge data sets
# 2. Random splits: Permutation test to construct empirical distribution
rule split:
    input:
        config["datasets"],
    output:
        merged = "merged/M0.rds",
        split1 = expand("splits/D1_s{num_splits}.rds", num_splits=range(config["num_splits"])),
        split2 = expand("splits/D2_s{num_splits}.rds", num_splits=range(config["num_splits"]))
    params:
        B = config["num_splits"],
        R = config["R"],
    log:
        "logs/split.log"
    shell:
        "{params.R} scripts/split.R -B {params.B} --log {log} --output {output.merged} {input}"

# 3. Run H-CBN2 (splits)
rule run_HCBN2_splits:
    input:
        "splits/D{file}_s{sample}.rds"
    output:
        "output/asa_D{file}_s{sample}.rds"
    params:
        sampling_mode = config["run_HCBN2"]["sampling_mode"],
        num_samples = config["run_HCBN2"]["num_samples"],
        iter_MCEM = config["run_HCBN2"]["iter_MCEM"],
        T0 = config["run_HCBN2"]["T0"],
        adap_rate = config["run_HCBN2"]["adap_rate"],
        iter_ASA = config["run_HCBN2"]["iter_ASA"],
        adaptive = "--asa" if config["run_HCBN2"]["adaptive"] else "", 
        seed = config["run_HCBN2"]["seed"],
        out_prefix = "output/D{file}_s{sample}_",
        R = config["R"],
    log:
        "logs/run_HCBN2_D{file}_s{sample}.log"
    threads:
        config["run_HCBN2"]["threads"]
    shell:
        "{params.R} scripts/run_HCBN2.R -L {params.num_samples} -s {params.sampling_mode} -i {params.iter_MCEM} -t {params.T0} -a {params.adap_rate} --iter_ASA {params.iter_ASA} {params.adaptive} --seed {params.seed} --out_prefix {params.out_prefix} --thrds {threads} --log {log} --output {output} {input}"

# 4. Run H-CBN2 (MLE)
rule run_HCBN2:
    input:
        get_input
    output:
        "output/asa_dataset{file}.rds"
    params:
        sampling_mode = config["run_HCBN2"]["sampling_mode"],
        num_samples = config["run_HCBN2"]["num_samples"],
        iter_MCEM = config["run_HCBN2"]["iter_MCEM"],
        T0 = config["run_HCBN2"]["T0"],
        adap_rate = config["run_HCBN2"]["adap_rate"],
        iter_ASA = config["run_HCBN2"]["iter_ASA"],
        adaptive = "--asa" if config["run_HCBN2"]["adaptive"] else "", 
        seed = config["run_HCBN2"]["seed"],
        out_prefix = "output/D{file}_",
        R = config["R"],
    log:
        "logs/run_HCBN2_D{file}.log"
    threads:
        config["run_HCBN2"]["threads"]
    shell:
        "{params.R} scripts/run_HCBN2.R -L {params.num_samples} -s {params.sampling_mode} -i {params.iter_MCEM} -t {params.T0} -a {params.adap_rate} --iter_ASA {params.iter_ASA} {params.adaptive} --seed {params.seed} --out_prefix {params.out_prefix} --thrds {threads} --log {log} --output {output} {input}"

# 5. Distance based tests
rule distance_based_test:
    input:
        D1 = "output/asa_dataset1.rds",
        D2 = "output/asa_dataset2.rds",
        D1_splits = expand("output/asa_D1_s{sample}.rds", sample=range(config["num_splits"])),
        D2_splits = expand("output/asa_D2_s{sample}.rds", sample=range(config["num_splits"])),
    output:
        "empirical_distribution.rds"
    params:
        B = config["num_splits"],
        R = config["R"],
    log:
        "logs/distance_based_test.log"
    shell:
        "{params.R} scripts/distance_based.R -B {params.B} --log {log} --output {output} {input}"


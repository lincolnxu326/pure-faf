# import basic packages
import pandas as pd
import sys
from snakemake.utils import validate
from cyvcf2 import VCF


# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"tumour_name": str})
)
if not samples["tumour_name"].is_unique:
    dups = samples["tumour_name"][samples["tumour_name"].duplicated()].unique()
    print("\nERROR: tumour_name values must be unique!")
    print("Duplicated tumour_name entries:", list(dups))
    sys.exit(1)

samples = samples.set_index("tumour_name", drop=False)
#exclude samples in the config exclude list
if "exclude_samples" in config:
    samples = samples[~samples["tumour_name"].isin(config["exclude_samples"])]
    
config["pval_threshold_high_AD"]=float(config["pval_threshold_high_AD"])
config["pval_threshold_low_AD"]=float(config["pval_threshold_low_AD"])

# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yaml")
validate(config, schema="../../config/schemas/config.schema.yaml")

# set config parameters
REF_GEN= config["reference_genome_fasta"]
ROQ= config["ROQ"]
AD= config["AD"]
DP= config["DP"]
VAF= config["VAF"]
PON_ALPHA= config["PON_alpha"]
AD_TIER=config["AD_tier_threshold"]
MSEC_PVAL_LOW=config["pval_threshold_low_AD"]
MSEC_PVAL_HIGH=config["pval_threshold_high_AD"]
ODIR= "results/" + config["odir"]


#set parameters based on reference genome version
if config['ref_genome_version'].lower()=='grch37': 
    FFPE_PON= "resources/FFPE_PON_hg19.bed"
    FFPE_PON_chr= "resources/FFPE_PON_hg19_chr.bed",
    SR_BED= "resources/trf_simple_repeat_hg19.bed"
    SR_BED_chr= "resources/trf_simple_repeat_hg19_chr.bed"
elif config['ref_genome_version'].lower()=='grch38':
    FFPE_PON= "resources/FFPE_PON_liftover_hg38.bed"
    FFPE_PON_chr= "resources/FFPE_PON_liftover_hg38_chr.bed"
    SR_BED= "resources/trf_simple_repeat_hg38.bed"
    SR_BED_chr= "resources/trf_simple_repeat_hg38_chr.bed"


def check_samples_in_vcf(vcf_path, tumour_sample, normal_sample):
    """
    Check that both tumour and normal samples exist in the VCF header.

    Exits the program with code 1 if either is missing.
    Returns True if both are present.
    """
    vcf = VCF(vcf_path)
    samples = vcf.samples

    missing = []
    if tumour_sample not in samples:
        missing.append(f"tumour sample '{tumour_sample}'")
    if normal_sample not in samples:
        missing.append(f"normal sample '{normal_sample}'")

    if missing:
        print(f"\nERROR: {vcf_path}")
        print("  Samples NOT found in the VCF file provided:", ", ".join(missing))
        print("  Samples found in the VCF file:", ", ".join(samples))
        sys.exit(1)

    return True

# use the function to check all samples in the sample sheet

for idx, row in samples.iterrows():
    check_samples_in_vcf(
        vcf_path=row["Mutect2_vcf"],
        tumour_sample=row["tumour_name"],
        normal_sample=row["normal_name"],
    )

# function to check if VCF uses 'chr' prefix in contig names
def vcf_uses_chr_prefix(vcf_path):
    """
    Determine whether the VCF uses 'chr' chromosome naming,
    based only on the first contig in the header.
    """
    vcf = VCF(vcf_path)
    contigs = list(vcf.seqnames)

    if not contigs:
        raise ValueError(f"No contigs in header for: {vcf_path}")

    first_contig = contigs[0]
    return first_contig.startswith("chr")

# set correct bed files based on VCF contig names
samples["chr"] = samples["Mutect2_vcf"].apply(vcf_uses_chr_prefix)
samples["PON"] = samples["chr"].map(
    lambda has_chr: FFPE_PON_chr if has_chr else FFPE_PON
)
samples["SR"] = samples["chr"].map(
    lambda has_chr: SR_BED_chr if has_chr else SR_BED
)


# functions to get sample specific parameters
def get_normal_name(wildcards):
    return samples.at[wildcards.tum, "normal_name"]

def get_mutect2_vcf(wildcards):
    return samples.at[wildcards.tum, "Mutect2_vcf"]

def get_deepsomatic_vcf(wildcards):
    return samples.at[wildcards.tum, "DeepSomatic_vcf"]

def get_tumour_bam(wildcards):
    return samples.at[wildcards.tum, "tumour_bam_cram"]

def get_ref_gen_fa(wildcards):
    return samples.at[wildcards.tum, "reference_genome_fasta"]

def get_sequencing_read_length(wildcards):
    return samples.at[wildcards.tum, "sequencing_read_length"]

def get_sequencing_adapter_1(wildcards):
    return samples.at[wildcards.tum, "sequencing_adapter_1"]

def get_sequencing_adapter_2(wildcards):
    return samples.at[wildcards.tum, "sequencing_adapter_2"]

def get_chr(wildcards):
    chr_status = "chr" if samples.at[wildcards.tum, "chr"] else "no_chr"
    return chr_status

def get_PON(wildcards):
    return samples.at[wildcards.tum, "PON"]

def get_SR(wildcards):
    return samples.at[wildcards.tum, "SR"]

def bsgenome_prereq():
    if config['ref_genome_version'].lower()=='grch37':
        return ".Rlib/hg19_download_complete.txt"
    return []
def bsgenome_prereq():
    ref = str(config["ref_genome_version"]).lower()
    if ref == "grch37":
        return f"{ODIR}/.Rlib/hg19_download_complete.txt"
    elif ref == "grch38":
        return f"{ODIR}/.Rlib/hg38_download_complete.txt"
    else:
        raise ValueError("ref_genome_version must be grch37 or grch38")



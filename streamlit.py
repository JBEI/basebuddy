from typing import Dict, List
import streamlit as st
import pandas as pd
from io import StringIO
import re

import numpy

from dnachisel import *
from Bio import Seq, SeqIO
from Bio.Restriction.Restriction_Dictionary import rest_dict
from streamlit_searchbox import st_searchbox
from thefuzz import process
from functools import reduce


# Note that this, which is gitignored, is also excluded from gcloud builds.
# For that reason, and cleanliness, that DB has been copied to a data/ directory.
# COCOPUTS_DB_FNAME = '220916_codon_analysis/o537-Refseq_species.tsv'
COCOPUTS_DB_FNAME = "data/cocoput_table.tsv"


# @st.experimental_memo
@st.experimental_singleton
def get_cocoput_organism_series():
    # df = pd.read_csv('220916_codon_analysis/220926_genome_codons.tsv', sep='\t', index_col=False)
    df = pd.read_csv(COCOPUTS_DB_FNAME, sep="\t", index_col=False)
    # df = pd.read_csv('220916_codon_analysis/o537-genbank_species.tsv', sep='\t', index_col=False)
    return pd.Series(
        df.apply(lambda r: f"{r['Species']} (TaxID: {r['Taxid']})", axis=1)
        .unique()
    )


@st.experimental_singleton
def get_cocoput_organism_list():
    return get_cocoput_organism_series().tolist()


def search_organisms(searchterm) -> List[str]:
    # print(f'search term: {searchterm}', flush=True)
    matches = get_cocoput_organism_series()
    matches = matches[reduce((lambda a, b: a & b), [matches.str.lower().str.contains(t.lower()) for t in searchterm.split()])]
    # for term in [t.lower() for t in searchterm.split()]:
    #     matches = matches[matches.str.contains(term)]
    if matches.shape[0] > 100:
        matches = matches.iloc[:100]
    # def matches_filter(item):
    #     item_lower = item.lower()
    #     for term in terms:
    #         if not term.contains(item_lower):
    #             return False
    #     return True
    # matches = [t for t in get_cocoput_organism_list() if matches_filter(t)]
    # matches = process.extract(searchterm, get_cocoput_organism_list(), 50)
    # matches = get_close_matches(searchterm, get_cocoput_organism_list(), 50)
    # print(matches, flush=True)
    return matches.tolist() #[match[0] for match in matches]


def get_taxid_from_cocoput_name(cocoput_name):
    rematch = re.match(r".*\(TaxID: (\d+)\)", cocoput_name)
    assert rematch, f"Somehow the cocoput name was poorly formatted, {cocoput_name}"
    return int(rematch.groups()[0])


AA_TO_CODON: Dict[str, List[str]] = {
    "*": ["TAA", "TAG", "TGA"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "C": ["TGC", "TGT"],
    "D": ["GAC", "GAT"],
    "E": ["GAA", "GAG"],
    "F": ["TTC", "TTT"],
    "G": ["GGA", "GGC", "GGG", "GGT"],
    "H": ["CAC", "CAT"],
    "I": ["ATA", "ATC", "ATT"],
    "K": ["AAA", "AAG"],
    "L": ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
    "M": ["ATG"],
    "N": ["AAC", "AAT"],
    "P": ["CCA", "CCC", "CCG", "CCT"],
    "Q": ["CAA", "CAG"],
    "R": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"],
    "S": ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
    "T": ["ACA", "ACC", "ACG", "ACT"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "W": ["TGG"],
    "Y": ["TAC", "TAT"],
}


def convert_cocoputs_table_to_dnachisel(
    codon_table_counts: dict,
) -> Dict[str, Dict[str, float]]:
    new_codon_table: Dict[str, Dict[str, float]] = {}
    for aa in AA_TO_CODON:
        new_codon_table[aa] = {}
        codon_sum: int = sum([codon_table_counts[codon] for codon in AA_TO_CODON[aa]])
        for codon in AA_TO_CODON[aa]:
            new_codon_table[aa][codon] = round(codon_table_counts[codon] / codon_sum, 3)
    return new_codon_table


@st.experimental_memo
def get_codon_table_for_taxid(taxid):
    df = pd.read_csv(COCOPUTS_DB_FNAME, sep="\t", index_col=False)
    subset = df[df.Taxid == taxid]
    row = subset[subset["# CDS"] == subset["# CDS"].max()].iloc[0]
    codons = [a + b + c for a in "ATCG" for b in "ATCG" for c in "ATCG"]
    codon_table = {k: v for k, v in row.to_dict().items() if k in codons}
    return convert_cocoputs_table_to_dnachisel(codon_table)


st.set_page_config(page_title="BaseBuddy")
# st.header('''BaseBuddy''')

c1, c2, c3 = st.columns((1, 5, 1))
with c2:
    st.image("resources/logo/png/logo-no-background.png")  # , width=500

st.write(
    """
An open-source based codon optimization tool for true transparency and reproducible results. 
"""
)


footer = """
<style>
  a:link , a:visited {
    color: blue;
    background-color: transparent;
    text-decoration: underline;
  }

  a:hover,  a:active {
    color: red;
    background-color: transparent;
    text-decoration: underline;
  }

  .footer {
    position: fixed;
    left: 0;
    bottom: 0;
    width: 100%;
    background-color: white;
    color: black;
    text-align: center;
  }
</style>
<div class="footer">
  <i>Powered by <a href="https://edinburgh-genome-foundry.github.io/DnaChisel/">DNA Chisel</a> and <a href="https://dnahive.fda.gov/dna.cgi?cmd=cuts_main">CoCoPUTs</a>.</i>
</div>
"""
st.markdown(footer, unsafe_allow_html=True)

col1, col2 = st.columns(2)

with col1:
    optimization_method = st.radio(
        "Optimization Method",
        ["use_best_codon", "match_codon_usage", "harmonize_rca"],
        key="visibility",
        # help='Choose the optimization method, harmonize_rca is recommended.'
        # label_visibility=st.session_state.visibility,
        # disabled=st.session_state.disabled,
        # horizontal=st.session_state.horizontal,
    )

with col2:
    st.caption('target organism')
    # target_organism = st.selectbox("target organism", get_cocoput_organism_list())
    target_organism = (st_searchbox(
        search_organisms,
        key="target_searchbox",
    ) or "")
    # print(f'TARGET: {target_organism}', flush=True)

    target_taxid = get_taxid_from_cocoput_name(target_organism) if target_organism else None
    target_coding_table = get_codon_table_for_taxid(target_taxid) if target_taxid else None

    if optimization_method == "harmonize_rca":
        source_organism = st.selectbox("source organism", get_cocoput_organism_list())
        source_taxid = get_taxid_from_cocoput_name(source_organism)
        source_coding_table = get_codon_table_for_taxid(source_taxid)
    else:
        source_coding_table = None


original_fasta_str = st.text_area(
    "Insert your native nucleotide sequence(s) in FASTA format:",
    ">example_sequence\nATGAGTAGT",
    help=f"It is important to only use the native coding sequence because some optimization methods calculate the relative codon adaptation (RCA) based on the inserted sequence.",
)
with StringIO(original_fasta_str) as fio:
    records = list(SeqIO.parse(fio, "fasta"))
if len(records) == 0:
    st.warning(f"Found zero valid records in text box, should have one.")
    st.stop()


# A list of constraints to impose, including cut site constraints.
cut_site_constraints = []

with st.expander("Advanced Settings"):
    homopolycols = st.columns(4)
    with homopolycols[0]:
        poly_a_maxlength = st.number_input(
            "Avoid poly As",
            value=9,
            min_value=1,
            max_value=15,
            step=1,
            help="Set the allowed maximum length of consecutive As.",
        )
    with homopolycols[1]:
        poly_t_maxlength = st.number_input(
            "Avoid poly Ts",
            value=9,
            min_value=1,
            max_value=15,
            step=1,
            help="Set the allowed maximum length of consecutive Ts.",
        )
    with homopolycols[2]:
        poly_c_maxlength = st.number_input(
            "Avoid poly Cs",
            value=6,
            min_value=1,
            max_value=15,
            step=1,
            help="Set the allowed maximum length of consecutive Cs.",
        )
    with homopolycols[3]:
        poly_g_maxlength = st.number_input(
            "Avoid poly Gs",
            value=6,
            min_value=1,
            max_value=15,
            step=1,
            help="Set the allowed maximum length of consecutive Gs.",
        )

    hairpin_c1, hairpin_c2 = st.columns((1, 1))
    with hairpin_c1:
        hairpin_stem_size = st.number_input(
            "Hairpin Stem Size",
            value=10,
            min_value=1,
            max_value=100,
            step=1,
            help="Set the allowed maximum hairpin stem size.",
        )
    with hairpin_c2:
        hairpin_window = st.number_input(
            "Hairpin Window",
            value=100,
            min_value=50,
            max_value=500,
            step=1,
            help="Define the minimum distance between hairpins.",
        )
    
    enforce_gc = st.columns(3)
    with enforce_gc[0]:
        gc_minimum = st.number_input(
            "Min. GC-content",
            value=0.3,
            min_value=0.1,
            max_value=0.9,
            step=0.01,
            help="Define the minium GC-content in a given window."
        )    

    with enforce_gc[1]:    
        gc_maximum = st.number_input(
            "Max. GC-content",
            value=0.75,
            min_value=0.1,
            max_value=0.9,
            step=0.01,
            help="Define the maximum GC-content in a given window."
        )

    with enforce_gc[2]:
        gc_window = st.number_input(
            "GC-content Window",
            value=50,
            min_value=10,
            max_value=200,
            step=1,
            help="Define the window to enforce the set GC-content."
        )
        

    restriction_sites = st.multiselect("Avoid cut sites",options=rest_dict,
        default=['BamHI', 'NdeI', 'XhoI', 'SpeI', 'BsaI']
        )
    for restricition_site in restriction_sites:
        cut_site_constraints.append(AvoidPattern(restricition_site+"_site")) 

    uniquify_kmers = st.columns(2)
    with uniquify_kmers[0]:
        kmers_value = st.number_input(
            "Uniquify All K-mers",
            value=10,
            min_value=1,
            max_value=20,
            step=1,
            help="The default value will ensure lower DNA synthesis diffculty due to less repetitive sequences. Changing to a lower value not recommended."
        )
    with uniquify_kmers[1]:
        kmers_rev = st.checkbox(
            "Include reverse complement",
            value=True,
            help="Uniquify the k-mers of the reverse complement."
        )
       
    randomize_numpy_seed = st.checkbox(
        "Numpy random generator",
        value=False,
        help="This ensures reproducible results. With this checked, an identical input seqeunce will always result in an identical output. However, unchecking this should not affect codon optimization result."
        )
    
        

    
    
    


# Do some input validation.
if not target_coding_table:
    st.warning("Must specify a target organism.")
    st.stop()
if optimization_method == "harmonize_rca":
    if not source_coding_table:
        st.warning("Must specify a source organism if using harmonize_rca method.")
        st.stop()
    codon_optimize_kwargs = {"original_codon_usage_table": source_coding_table}
else:
    codon_optimize_kwargs = {}

constraints_logs = []
objectives_logs = []
recodings = []
try:
    for record in records:

        if randomize_numpy_seed:
            numpy.random.seed(None)
        else:
            numpy.random.seed(123)
        problem = DnaOptimizationProblem(
            sequence=record,
            constraints=[
                UniquifyAllKmers(kmers_value, include_reverse_complement=kmers_rev),
                AvoidHairpins(
                    stem_size=hairpin_stem_size, hairpin_window=hairpin_window
                ),
                AvoidPattern(str(poly_a_maxlength) + "xA"),
                AvoidPattern(str(poly_t_maxlength) + "xT"),
                AvoidPattern(str(poly_c_maxlength) + "xC"),
                AvoidPattern(str(poly_g_maxlength) + "xG"),
                EnforceGCContent(mini=gc_minimum, maxi=gc_maximum, window=gc_window),
                EnforceTranslation(),
            ] + cut_site_constraints,
            objectives=[
                CodonOptimize(
                    method=optimization_method,
                    codon_usage_table=target_coding_table,
                    **codon_optimize_kwargs,
                )
            ],
        )

        problem.max_random_iters = 10000
        problem.resolve_constraints()
        problem.optimize()

        constraints_logs.append(problem.constraints_text_summary())
        objectives_logs.append(problem.objectives_text_summary())
        recodings.append(problem.sequence)

except Exception as e:
    st.warning(e)
    st.stop()


with st.expander("DNA Chisel Logs"):
    for log in constraints_logs + objectives_logs:
        st.text(log)

for record, recoding in zip(records, recodings):
    if optimization_method == 'harmonize_rca':
        notes = f'method: {optimization_method}, source_taxid: {source_taxid}, target_taxid: {target_taxid}'  
    else:
        notes = f'method: {optimization_method}, target_taxid: {target_taxid}'  
    st.text(f">{record.id} ({notes})\n{recoding}")

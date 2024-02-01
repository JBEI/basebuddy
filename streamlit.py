from typing import Dict, List
import streamlit as st
import pandas as pd
from io import StringIO
import re
from datetime import datetime


import numpy

from dnachisel import *
from Bio import Seq, SeqIO
from Bio.Restriction.Restriction_Dictionary import rest_dict
from streamlit_searchbox import st_searchbox
from thefuzz import process
from functools import reduce
from streamlit_modal import Modal
import streamlit.components.v1 as components
from st_keyup import st_keyup

# Note that this, which is gitignored, is also excluded from gcloud builds.
# For that reason, and cleanliness, that DB has been copied to a data/ directory.
# COCOPUTS_DB_FNAME = '220916_codon_analysis/o537-Refseq_species.tsv'
COCOPUTS_DB_FNAME = "data/cocoput_table.tsv"
COCOPUTS_INDEX_FNAME = "data/cocoput_index.csv"


# @st.experimental_memo
@st.experimental_singleton
def get_cocoput_organism_index():
    # df = pd.read_csv('220916_codon_analysis/220926_genome_codons.tsv', sep='\t', index_col=False)
    df = pd.read_csv(COCOPUTS_INDEX_FNAME, index_col=False)
    # df = pd.read_csv('220916_codon_analysis/o537-genbank_species.tsv', sep='\t', index_col=False)
    return pd.Series(
        df.apply(lambda r: f"{r['Species']} (TaxID: {r['Taxid']})", axis=1)
        .unique()
    )


@st.experimental_singleton
def get_cocoput_organism_list():
    return get_cocoput_organism_index().tolist()


def search_organisms(searchterm) -> List[str]:
    # print(f'search term: {searchterm}', flush=True)
    print(f'Searching for {searchterm}', flush=True)
    matches = get_cocoput_organism_index()
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

st.write("""
**Instructions:**

**1.) Choose your desired optimization method.**

**2.) Specify your target organism.**

**3.) Paste your coding sequence(s) into the textbox below and press Enter or upload a (multi-)FASTA file.**

**4.) Click the "Download Result(s)" button to generate a FASTA file with your optimized sequence(s).**

Check out the advanced settings if you have special requirements for your optimized sequence 
(e.g. avoiding certain cut sites or patterns, or lowering synthesis difficulty). 

For further questions or concerns please visit https://github.com/JBEI/basebuddy
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
        "**Optimization Method**",
        ["use_best_codon", "match_codon_usage", "harmonize_rca"],
        key="visibility",
        help='Choose a codon optimization method. See details on the [DNAChisel website](https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html#codon-optimization-specifications).'
        # label_visibility=st.session_state.visibility,
        # disabled=st.session_state.disabled,
        # horizontal=st.session_state.horizontal,
    )
    database = st.radio(
        "**Codon Usage Database**",
        ["CoCoPUTs", "Kazusa"],
        key="database",
        help='Choose a codon usage table database. CoCoPUTs is considered more accurate and up-to-date.'
        # label_visibility=st.session_state.visibility,
        # disabled=st.session_state.disabled,
        # horizontal=st.session_state.horizontal,
    )

with col2:
    # target_searchterm = st_keyup('Search for a target organism:')
    # st.dataframe(
    #     pd.DataFrame({'organism': search_organisms(target_searchterm)}, index=None),
    #     use_container_width=True
    # )
    # hide_table_row_index = """
    # <style>
    # tr:first-child {display:none}
    # tbody th {display:none}
    # </style>
    # """
    # # Inject CSS with Markdown
    # st.markdown(hide_table_row_index, unsafe_allow_html=True)
    # # st.table(search_organisms(target_searchterm))
    # df = pd.DataFrame(search_organisms(target_searchterm))
    # st.dataframe(df)

    # styler = df.style.hide_index()
    # st.write(styler.to_html(), unsafe_allow_html=True)
    # st.table(search_organisms(target_searchterm))

    # target_organism = st.selectbox("**Target Organism**", get_cocoput_organism_list())
    st.caption('**Target Organism**')
    target_organism = (st_searchbox(
        search_organisms,
        key="target_searchbox",
    ) or "")
    # print(f'TARGET: {target_organism}', flush=True)

    target_taxid = get_taxid_from_cocoput_name(target_organism) if target_organism else None
    target_coding_table = get_codon_table_for_taxid(target_taxid) if target_taxid else None

    if optimization_method == "harmonize_rca":
        st.caption("**Source Organism**")
        source_organism =  (st_searchbox(
            search_organisms,
            key="source_organism",
        ) or "")  #st.selectbox("**Source Organism**", get_cocoput_organism_list())
        # source_organism = st.selectbox("**Source Organism**", get_cocoput_organism_list())

        source_taxid = get_taxid_from_cocoput_name(source_organism) if source_organism else None
        source_coding_table = get_codon_table_for_taxid(source_taxid) if source_taxid else None
    else:
        source_coding_table = None


paste_tab, upload_tab = st.tabs(['Paste', 'FASTA Upload'])

with paste_tab:
    original_fasta_str = st.text_area(
        "**Insert your native nucleotide sequence(s) in FASTA format:**",
        ">example_sequence1\nATGAGTAGT\n>example_sequence2\nATGGTGAATTTG",
        help=f"It is important to only use the native coding sequence because some optimization methods calculate the relative codon adaptation (RCA) based on the inserted sequence.",
    )
with upload_tab:
    uploaded_fasta_file = st.file_uploader(
        "**Upload your (multi-)FASTA file here:**", type=[".fasta",".fa"]
    )
    if uploaded_fasta_file is not None:
        text_io = uploaded_fasta_file.read().decode("UTF-8")
        original_fasta_str = text_io
    
  

with StringIO(original_fasta_str) as fio:
    records = list(SeqIO.parse(fio, "fasta"))
if len(records) == 0:
    st.warning(f"Found zero valid records in textbox, should have one.")
    st.stop()

def process_substring(substring):
    # Convert substring to uppercase, strip whitespace
    return substring.strip().upper()

def validate_input_string(input_string):
    allowed_chars = {'A', 'T', 'G', 'C', ','}
    for char in input_string:
        if char.upper() not in allowed_chars and char != ',' and not char.isspace():
            return False
    return True

# A list of constraints to impose, including cut site constraints.
cut_site_constraints = []
custom_pattern_constraints = []

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
        

    restriction_sites = st.multiselect(
        "Avoid cut sites",
        options=rest_dict,
        default=['BamHI', 'NdeI', 'XhoI', 'SpeI', 'BsaI']
        )
    for restricition_site in restriction_sites:
        cut_site_constraints.append(AvoidPattern(restricition_site+"_site")) 

    # Get Avoid custom pattern string from the user
    input_string = st.text_input(
        "Avoid custom pattern(s)",
        placeholder="e.g. atta, gggtttaaa, ...",
        help='''Adding a single DNA sequence or multiple sequences divided by commas will avoid that specific sequence pattern on both strands. Only As, Ts, Gs, Cs, or commas are allowed (not case-sensitive). For example, the input:"atta, agggt" will avoid the sequence pattern "ATTA" and "AGGGT" in the optimized output.'''
    )

    
    if input_string.strip():  # Check if input_string is not empty after stripping whitespace
        if validate_input_string(input_string):
            # Split the input string by commas
            input_sequences = input_string.split(',')
            
            # Process each input sequence and create constraints
            for sequence in input_sequences:
                processed_sequence = process_substring(sequence)
                if processed_sequence:  # Check if the processed sequence is non-empty
                    custom_pattern_constraints.append(AvoidPattern(processed_sequence))  # Assuming AvoidPattern accepts the processed sequence
        else:
            st.warning("Error: Input string contains invalid characters. Only 'As', 'Ts', 'Gs', 'Cs', or commas are allowed.")
    



    uniquify_kmers = st.columns(2)
    with uniquify_kmers[0]:
        kmers_value = st.number_input(
            "Uniquify All K-mers",
            value=10,
            min_value=1,
            max_value=20,
            step=1,
            help="This will ensure lower DNA synthesis diffculty due to less repetitive codons. Changing to a lower value might further lower synthesis difficulty but can also result in an unsolvable constraint."
        )
    with uniquify_kmers[1]:
        kmers_rev = st.checkbox(
            "Include reverse complement",
            value=True,
            help="Uniquify the k-mers of the reverse complement as well."
        )
       
    randomize_numpy_seed = st.checkbox(
        "Numpy random generator",
        value=False,
        help="This ensures reproducible results. With this unchecked, an identical input sequence will always result in an identical output. However, checking this should not affect codon optimization result."
        )
    
        

    
    
    


# Do some input validation.
if not target_coding_table:
    st.warning("Please specify a target organism.")
    st.stop()


if database == 'CoCoPUTs':
    codon_optimize_kwargs = {
        "codon_usage_table": target_coding_table
    }
    if optimization_method == "harmonize_rca":
        if not source_coding_table:
            st.warning("Please specify a source organism if using the harmonize_rca method.")
            st.stop()
        codon_optimize_kwargs["original_codon_usage_table"] = source_coding_table

elif database == 'Kazusa':
    codon_optimize_kwargs = {
        "species": target_taxid
    }
    if optimization_method == "harmonize_rca":
        if not source_taxid:
            st.warning("Please specify a source organism if using the harmonize_rca method.")
            st.stop()
        codon_optimize_kwargs["original_species"] = source_taxid
        
else:
    raise KeyError(f'Unrecognized database {database}. Maybe the database dropdown was reformatted?')

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
                
            ] + cut_site_constraints + custom_pattern_constraints,
            objectives=[
                CodonOptimize(
                    method=optimization_method,
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

result_list=[]
for record, recoding in zip(records, recodings):
    if optimization_method == 'harmonize_rca':
        notes = f'method: {optimization_method}, source_taxid: {source_taxid}, target_taxid: {target_taxid}'  
    else:
        notes = f'method: {optimization_method}, target_taxid: {target_taxid}'
    
    record_result = f">{record.id} ({notes}) Please cite us :) \n{recoding}\n"
    result_list.append(record_result)


st.download_button(label="**Download result(s)**", 
                    data = "".join(str(j) for j in result_list),
                    file_name = datetime.now().strftime("%Y%m%d-%I%M%S%p_") + "BaseBuddy_results" + ".fasta",
                    )   
st.text_area(
    "Or copy your result(s) from this textbox:","".join(str(j) for j in result_list)
)

with st.expander("Optimization Logs"):
    for log in constraints_logs + objectives_logs:
        st.text(log)
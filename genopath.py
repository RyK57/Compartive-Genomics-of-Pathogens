import streamlit as st
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
from Bio import pairwise2


# Set Streamlit page configuration
st.set_page_config(page_title="Comparative Genomic Analysis")

# Primary accent for interactive elements
primaryColor = '#7792E3'

# Background color for the main content area
backgroundColor = '#273346'

# Background color for sidebar and most interactive widgets
secondaryBackgroundColor = '#B9F1C0'

# Color used for almost all text
textColor = '#FFFFFF'

# Font family for all text in the app, except code blocks
font = "sans serif"

def read_fasta(file):
    # Reading the contents of the file
    contents = file.getvalue().decode("utf-8")
    # Splitting the contents into lines
    lines = contents.split("\n")
    # Removing any empty lines
    lines = [line for line in lines if line.strip()]
    # Parsing the FASTA sequences
    sequences = {}
    max_length = 0
    for line in lines:
        if line.startswith(">"):
            current_sequence = line[1:]
            sequences[current_sequence] = ""
        else:
            sequences[current_sequence] += line

            # Keeping track of the maximum sequence length
            if len(sequences[current_sequence]) > max_length:
                max_length = len(sequences[current_sequence])

    # Adding padding to ensure all sequences are the same length
    for sequence_name in sequences:
        sequence = sequences[sequence_name]
        padding = max_length - len(sequence)
        sequences[sequence_name] = sequence + "-" * padding

    # Converting the sequences dictionary to a list of SeqRecord objects
    seq_records = [SeqRecord(Seq(seq), id=name) for name, seq in sequences.items()]
    # Aligning the sequences and return the result
    return align_sequences(seq_records)


def align_sequences(seq_records):
    # Using the first sequence as the template for the alignment
    template_seq = seq_records[0].seq
    # Creating a list of sequences to be aligned
    seqs_to_align = [str(seq_record.seq) for seq_record in seq_records]
    # Aligning the sequences
    alignment = MultipleSeqAlignment(seq_records)
    # Returning the alignment
    return alignment

def build_comparision_graph(alignment):
    # Creating a distance calculator object
    calculator = DistanceCalculator('identity')
    # Calculating the distances between sequences in the alignment
    distance_matrix = calculator.get_distance(alignment)
    # Creating a distance-based tree constructor object
    constructor = DistanceTreeConstructor()
    # Building the tree
    tree = constructor.upgma(distance_matrix)
    return tree


# Sidebar
st.sidebar.title("Comparative Genomic Analysis")
st.sidebar.write("-Rithvik Sabnekar")
st.markdown("NCBI Virus Genome Browser - https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239")
file_url = 'https://drive.google.com/file/d/1vwHW4NcsU6tfJPbKVr1CGiZLtz5NpgtY/view?usp=sharing'

# link to download the example sequence file
st.markdown(f'<a href="https://drive.google.com/file/d/1coCSpNrDI599WICcsCj5cs6C8U9Hb46g/view?usp=sharing" download>Example Input file (Carpodacus Mexicanus)</a>', unsafe_allow_html=True)


# File Upload
uploaded_file = st.sidebar.file_uploader("Upload a FASTA file", type=["txt"])

if uploaded_file is not None:
    # Read and display the uploaded file
    alignment = read_fasta(uploaded_file)

    # Alignment
    st.header("Alignment done")

    # Consensus Sequence
    st.subheader("Consensus Sequence derived")
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()

    # Distance Matrix
    st.subheader("Distance Matrix calculated")
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    st.write(distance_matrix)

    # Comparative Genomic Analysis
    st.subheader("Comparative Genomic Analysis")
    tree = build_comparision_graph(alignment)

    # Plot the tree
    plt.rcParams["figure.figsize"] = [10, 7]  # Customize the plot size
    Phylo.draw(tree)
    plt.title("Comparative Genomic Analysis")
    plt.xlabel("Distance")
    plt.ylabel("Sequences")
    st.pyplot(plt)

    # Sequence Length Distribution
    st.subheader("Sequence Length Distribution")
    sequence_lengths = [len(record.seq) for record in alignment]
    plt.hist(sequence_lengths, bins=20)
    plt.title("Sequence Length Distribution")
    plt.xlabel("Sequence Length")
    plt.ylabel("Count")
    st.pyplot(plt)

    # Pairwise Sequence Similarity Heatmap
    st.subheader("Pairwise Sequence Similarity Heatmap")
    similarity_matrix = []
    for i, record1 in enumerate(alignment):
        similarity_row = []
        for j, record2 in enumerate(alignment):
            alignment_score = pairwise2.align.globalxx(record1.seq, record2.seq, one_alignment_only=True)[0].score
            similarity_row.append(alignment_score)
        similarity_matrix.append(similarity_row)

    plt.imshow(similarity_matrix, cmap="hot", interpolation="nearest")
    plt.colorbar(label="Similarity Score")
    plt.title("Pairwise Sequence Similarity Heatmap")
    plt.xlabel("Sequence")
    plt.ylabel("Sequence")
    st.pyplot(plt)

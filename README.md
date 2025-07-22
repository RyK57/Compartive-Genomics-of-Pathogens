# Compartive-Genomics-of-Pathogens
Streamlit app to visualize relation and characterisitics between various viruses using their genomic data

Function:
The application allows users to upload a FASTA file containing genomic sequences of pathogens and provides several analysis outputs:

Alignment: The uploaded sequences are aligned using the Biopython library's MultipleSeqAlignment and align_sequences functions.
Consensus Sequence: It calculates the consensus sequence derived from the alignment using AlignInfo.SummaryInfo and dumb_consensus function.
Distance Matrix: The code calculates the distance matrix between sequences in the alignment using the DistanceCalculator class from Biopython.
Comparative Genomic Analysis: It builds a distance-based tree using the DistanceTreeConstructor and upgma function from the Biopython Phylo.TreeConstruction module. The tree represents the relationships between the different sequences.
Plotting: The application visualizes the comparative genomic analysis by plotting the constructed tree, sequence length distribution, and pairwise sequence similarity heatmap using the matplotlib.pyplot library.
Technical details of the code:

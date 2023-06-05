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

The application is built using the Streamlit framework, which allows for interactive web applications in Python.
It utilizes various modules and functions from the Biopython library, such as Bio.Align, Bio.SeqRecord, Bio.Seq, Bio.Phylo.TreeConstruction, and Bio.pairwise2, for sequence alignment, consensus calculation, distance calculation, tree construction, and pairwise sequence similarity.
The user interface is designed using Streamlit's functions, such as st.sidebar, st.markdown, and st.file_uploader, to create a sidebar, display text and links, and upload the FASTA file.
The analysis results are displayed using Streamlit's st.header, st.subheader, st.write, and st.pyplot functions.
Importance in performing comparative genomics:
Comparative genomic analysis is crucial for understanding the relationships, evolutionary patterns, and functional implications of various pathogens. By aligning and analyzing genomic sequences, researchers can identify conserved regions, discover genetic variations, infer evolutionary relationships, and gain insights into the structure and function of pathogens. The provided application simplifies the process by providing a user-friendly interface for performing these analyses and visualizing the results, such as the constructed tree, sequence length distribution, and pairwise sequence similarity heatmap. Researchers can use this application to gain a better understanding of pathogen genomics, identify potential virulence factors, track the spread of pathogens, and make informed decisions regarding treatments and preventive measures.

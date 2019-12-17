This is a collection of some of my solutions to problems posed in the course Bioinformatics Algorithms and Models (Professor David Bernick, UCSC, 2019). My algorithms are implemented in python, many of them are formatted for bioinformatics learning platform Rosalind. Detailed docstrings are included in each program.

Algorithms covered in these programs include:

- Stochastic Processes based on Shannon’s Entropy for Promotor Sequence Discovery
- Eulerian Path from Directed Graph
- DeBruijn Graphing
- String (quasi-genome) Assembly
- Topological Sorting (Directed Acyclic Graphs)
- Viterbi Learning for HMM Path Maximization.
- Branch and Bound Algorithm for Cyclic Peptide Reconstruction


Executables:

python3 randomizedMotifSearch.py -i int -p int -k int -r (opt) -g (opt) -m (opt) < ./input_files/pcaCrisprs

python3 missingMotif.py -minMotif 3 -maxMotif 8 -cutoff 0.0 < ./input_files/Zm4-genomic.fna 

python3 deBruijnGraphFromString.py  < ./input_files/deBruijn.txt

python3 stringReconstruction.py < ./input_files/stringReconstruction_in.txt

python3 longestPathInDAG.py < ./input_files/longestPathInDAG_in.txt

python3 viterbiLearning.py < ./input_files/viterbiLearning_in.txt

python3 cyclicPeptideReconstruction.py < ./input_files/experimentalSpectrum.txt
# Conservation-analysis
A pipeline for analyzing amino acid conservation

## Introduction
Analyzing amino acid conservation can help identify key residues or regions of interest in your given protein. This pipeline provides a quick way of assessing the conservation of each residue in your protein of interest.

The script uses 8 threads and 4GB RAM.

## Installation
Download the github repository
```
git clone https://github.com/nickatirwin/Conservation-analysis.git
```
If conda is not installed, install it
```
# got to home directory
cd ~
# download installer
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
# install conda
bash Miniforge3-$(uname)-$(uname -m).sh
```

Install dependencies using conda
```
conda env create -f environment.yml
# or
mamba env create -f environment.yml
```
Activate the conda environment
```
conda activate conservation
```
## Usage
The script uses pre-compiled datasets available on the CLIP cluster. These include a pan-eukaryotic database (eukaryota) and a archaeplastida database (including plants, green algae, red algae, and glaucophytes).

Create a fasta file for your query protein
```
# eg. AtDDM1.fasta
>AtDDM1
MVSLRSRKVIPASEMVSDGKTEKDASGDSPTSVLNEEENCEEKSVTVVEEEILLAKNGDSSLISEAMAQEEEQLLKLREDEEKANNAGSAVAPNLNETQFTKLDELLTQTQLYSEFLLEKMEDITINGIESESQKAEPEKTGRGRKRKAASQYNNTKAKRAVAAMISRSKEDGETINSDLTEEETVIKLQNELCPLLTGGQLKSYQLKGVKWLISLWQNGLNGILADQMGLGKTIQTIGFLSHLKGNGLDGPYLVIAPLSTLSNWFNEIARFTPSINAIIYHGDKNQRDELRRKHMPKTVGPKFPIVITSYEVAMNDAKRILRHYPWKYVVIDEGHRLKNHKCKLLRELKHLKMDNKLLLTGTPLQNNLSELWSLLNFILPDIFTSHDEFESWFDFSEKNKNEATKEEEEKRRAQVVSKLHGILRPFILRRMKCDVELSLPRKKEIIMYATMTDHQKKFQEHLVNNTLEAHLGENAIRGQGWKGKLNNLVIQLRKNCNHPDLLQGQIDGSYLYPPVEEIVGQCGKFRLLERLLVRLFANNHKVLIFSQWTKLLDIMDYYFSEKGFEVCRIDGSVKLDERRRQIKDFSDEKSSCSIFLLSTRAGGLGINLTAADTCILYDSDWNPQMDLQAMDRCHRIGQTKPVHVYRLSTAQSIETRVLKRAYSKLKLEHVVIGQGQFHQERAKSSTPLEEEDILALLKEDETAEDKLIQTDISDADLDRLLDRSDLTITAPGETQAAEAFPVKGPGWEVVLPSSGGMLSSLNS
```
Run the script
```
# choose whether to use the eukaryota or archaeplastida database (these are present within the CLIP cluster)
python conservation.py AtDDM1.fasta eukaryota
python conservation.py AtDDM1.fasta archaeplastida
```
Note: The database could also be customized.
```
# eg. custom.fasta
# make a diamond database
diamond makedb --in custom.fasta --db custom.dmnd
python conservation.py AtDDM1.fasta datasets/custom
```

## Output
Results will be output into a results directory. The script will output formatted plots showing the conservation score inferred from the phylogenetic analysis. Conservation is calculated as 1/rate of each site. A table compiling all results is created (ends with .tsv).

The alignment from which the conservation scores was calculated is in the results directory and called: query.select.fasta.aln.trim.

## Methods
The script takes the query protein as input and searches the designated database using DIAMOND BLASTP (e-value < 1e-5, query-coverage > 50, max hits: 1000). The hits are then aligned using MAFFT, trimmed based on the query sequence, and an initial phylogeny is inferred using IQTREE (model: LG+I+G4, fast mode). The tree is then loaded in ETE and the tree is traversed from the query protein. The clade with sequences closest to 200 will be selected, and the resulting sequences are extracted. The extracted sequences are then re-aligned using MAFFT, trimmed based on the query, and a second phylogeny is ran (using the same parameters, as well as -wsr for rate reporting). The resulting rates are plotted and reported alongside the Shannon entropy and sequence distance (the median distance between the query amino acid and the other sequences based on BLOSUM62).  




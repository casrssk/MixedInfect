# MixedInfect
#### A strain identification tool for *Chlamydia trachomatis* whole genomes sequencing data using k-mers containing SNP information.

#### For contacting the developer and issue reports please go to https://github.com/casrssk/MixedInfect/issues

## Getting Started
##### These instruction will guide you through the process of deploying MixedInfect on your local computer. Tested to work on MAC OSX and Ubuntu.

### Setup
##### To used MixedInfect the following dependencies have to be satisfied
- Python3.5 or higher
- biopython (v1.78)
- harvesttools (v1.2)
- mash (v2.2)
- parsnp (v1.5.3)
- pyahocorasick (v1.4.0)
- samtools (v1.11)

##### you can also create a virtual environment and install the dependencies using mixed_infect.yml

##### To install MixedInfect, run:
- git clone https://github.com/casrssk/MixedInfect.git
- cd MixedInfect
- if you are using conda make sure that you activate the environment.

### Usage

##### run preprocess.py to build the k-mer database
- preprocess.py -d {dir path for referenec genomes} --ref_genome {path for referenec genome}

##### run freq_mat.py to generate frequency matrix
- freq_mat.py --ref_genome {path for referenec genome}

##### run mixedinfect.py to get the strain diversity and their relative abundances
- mixedinfect.py --path {path for sample fastq files} --ref_genome {path for referenec genome}

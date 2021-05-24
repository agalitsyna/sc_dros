sc_dros is an illustrative code for the paper where we unravel the story on _Drosophila_ chromatin in single cells: 

   ```
   Ulianov, Sergey V., Vlada V. Zakharova, Aleksandra A. Galitsyna, Pavel I. Kos, Kirill E. Polovnikov, Ilya M. Flyamer, Elena A. Mikhaleva et al. "Order and stochasticity in the folding of individual Drosophila genomes." Nature Communications 12, no. 1 (2021): 1-17. DOI: [10.1038/s41467-020-20292-z](https://doi.org/10.1038/s41467-020-20292-z)
   ```

This code is focused on initial steps of data processing. 
   
By default, it is supposed to work for _Drosophila_ data from the paper (see the notes, however). It may be adapted for other datasets (I've already tested it on Flyamer et al. 2017 and Gassler et al. 2017 data). In case of questions, don't hesitate to contact the maintainer: Aleksandra Galitsyna (Aleksandra.Galitsyna at skoltech.ru)

   

#### Preparing the environment

We recommend setting up conda environment named ```sc_dros``` with all the requirements.

List of dependencies: 

- snHi-C and Hi-C data processing: 

    bwa
  
    cooler 
    
    pairtools, version https://github.com/agalitsyna/pairtools.git
   
    lavaburst    

- Pyhon 3.7, libraries: 

    pandas

    numpy 
    
    seaborn


#### Preparing the dataset

The files needed:

- Files with information about the genome: 

    - Reduced chromosomes sizes,
    
      Path: ```data/GENOME```
    
    - FASTA file with genome and bwa index,
    
        Path: ```data/GENOME```
    
        Run: bash script ```scripts/00_prepare_data/001_prepare_genome.sh```
    
    - FASTQ files with snHi-C data,
        
        PATH: ```data/FASTQ```
        
        Run python code: ```scripts/00_prepare_data/002_download_data.py```
        
        The result of this step should be a set of FASTQ files with names formatted as follows: 
        
      ```
      {Cell_name}_{Replicate}_R1.fastq.gz
      {Cell_name}_{Replicate}_R2.fastq.gz
      ```

        _Note 1._ I recommend the specialized version of GEOParse.
        Dependencies: GEOparse with parallel download option: https://github.com/agalitsyna/GEOparse.git
        To install GEOparse, run: 
        
        ```bash
        git clone https://github.com/agalitsyna/GEOparse.git
        cd GEOparse
        pip install -e .
        ```
        
        _Note 2._ Note 2. GEOParse for many GEO IDs but may fail for some SRA entries. This is because GEOParse downloads the data from SRA FTP, [decommissioned at the end of 2019](https://ncbiinsights.ncbi.nlm.nih.gov/2019/10/17/users-of-the-sra-ftp-site-try-the-sra-toolkit/), and files might not be available through GEOparse.
      
        _Note 3._ If GEOParse does not install in your environment, or SRA fails to be downloaded from ftp, you may use manual **download from GEO**.
        FASTQs can be found in [GEO entry GSE131811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131811).
        

#### Raw data processing

Raw data processing includes the steps from raw FASTQ files to processed COOL files for both snHi-C and bulk Hi-C.

- Reads mapping and BAM files to PAIR parse
    
```bash
bash scripts/01_data_mapping/010_run_mapping.sh
bash scripts/01_data_mapping/011_parse_population.sh
```

- PAIR files to COOL conversion
```bash 
bash 012_run_pairsam2cooler.sh
```

This should result in a set of COOL files: ```data/COOL/```

_Note._ If there are problems with running the scripts/accessig the data, Drosophila snHi-C COOLs can be alternatively found as supplementary files in [GEO entry GSE131811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131811).

#### TAD calling

TAD calling is performed with lavaburst package and includes the steps of scanning a wide range of gamma parameter and selection of a single set of TADs and sub-TADs. 

Working directory: ```scripts/02_tad_calling/```

Script example: 021_run_FindOptimalGamma.sh


### Remarks

This folder is maintained for demonstration of initial steps of snHi-C data processing.

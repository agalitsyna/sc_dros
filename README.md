
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
        _Note 1._ Sometimes SRA storage experiences internal problems with ftp, and files might not be available through GEOparse. 

        _Note 2._ I recommend the specialized version of GEOParse.
        Dependencies: GEOparse with parallel download option: https://github.com/agalitsyna/GEOparse.git
        To install GEOparse, run: 
        
        ```bash
        git clone https://github.com/agalitsyna/GEOparse.git
        cd GEOparse
        pip install -e .
        ```
      
        _Note 3._ If GEOParse does not install in your environment, or SRA fails to be downloaded from ftp, you may use the direct download from GEO.
        The resulting FASTQs can be found in [GEO entry GSE131811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131811).

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

_Note._ The resulting COOLs can be alternatively found as supplementary files in [GEO entry GSE131811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131811).
#### TAD calling

TAD calling is performed with lavaburst package and includes the steps of scanning a wide range of gamma parameter and selection of a single set of TADs and sub-TADs. 

Working directory: ```scripts/02_tad_calling/```

Script example: 021_run_FindOptimalGamma.sh


### Remarks

This folder is maintained for demonstration of initial steps of snHi-C data processing.

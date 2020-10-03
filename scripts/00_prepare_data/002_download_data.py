# Creates the folder ../../data/FASTQ/ from GEO ID GSE131811

from sys import argv
if len(argv)<2:
    print("Please, specify your e-mail for GEO download as the first parameter of this script.")
    exit()

import GEOparse

fastq_dump_options={'split-files': None,
    'read-filter': 'pass',
    'dumpbase': None,
    'gzip': None}

geo_id = 'GSE131811'

gse = GEOparse.get_GEO(geo=geo_id, destdir="../../data/FASTQ/TMP_SOFT")
gsms = gse.gsms
downloaded_paths = gse.download_SRA(argv[1], directory="../../data/FASTQ/",
        filetype='fastq', 
        fastq_dump_options=fastq_dump_options, 
        nproc=5, 
        silent=False)

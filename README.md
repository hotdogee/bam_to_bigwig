# bam_to_bigwig

Convert your alignment files in BAM format into coverage files in bigWig format

# Usage
```
bam_to_bigwig.py BAM_FILE [BAM_FILE ...]
[-o, --bigwig_filename=<output file name>
 -t, --tempfile
 -k, --keep-tempfile
 -s, --ignore-secondary
 -q, --ignore-qc-fail
 -d, --ignore-optical-pcr-duplicate
 -u, --ignore-supplementary]

--bigwig_filename: if not given will save to <BAM file prefix>.bigwig.
--tempfile: If this is set, will use tempfile library to generate file names instead of using 
            <BAM file prefix>.wig and <BAM file prefix>.sizes.
--keep-tempfile: If this is set, will not delete the intermediate files <BAM file prefix>.wig 
                 and <BAM file prefix>.sizes.
--ignore-secondary: If this is set, rsem-bam2wig will ignore alignments with the 
                    "secondary alignment" flag bit 0x100 set.
--ignore-qc-fail: If this is set, rsem-bam2wig will ignore alignments with the 
                  "not passing quality controls" flag bit 0x200 set.
--ignore-optical-pcr-duplicate: If this is set, rsem-bam2wig will ignore alignments with the 
                                "PCR or optical duplicate" flag bit 0x400 set.
--ignore-supplementary: If this is set, rsem-bam2wig will ignore alignments with the 
                        "supplementary alignment" flag bit 0x800 set.
```

# Dependencies
* Python 2.7 on CentOS 6 x86_64 (may work on windows if you can get the other dependencies to compile)
  * pysam (https://github.com/pysam-developers/pysam)
* wigToBigWig from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/)
* rsem-bam2wig from RSEM (http://deweylab.biostat.wisc.edu/rsem/)

# Installation
* pysam
  1. Install build dependencies
    * CentOS: yum install python-devel zlib-devel
    * Ubuntu: apt-get install python-dev zlib1g-dev
  2. Install using pip
    * pip install pysam
* wigToBigWig
  1. Download binary: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
  2. Make sure wigToBigWig is in a directory in your $PATH
    * cp wigToBigWig /usr/local/bin
* rsem-bam2wig
  1. Download the source code: http://deweylab.biostat.wisc.edu/rsem/
  2. Compile
    * make
  3. Make sure rsem-bam2wig is in a directory in your $PATH
    * cp rsem-bam2wig /usr/local/bin

# Authors
```
Han Lin, 2014.
Biomedical Electronics and Bioinformatics, National Taiwan University
National Agricultural Library, ARS, USDA
```

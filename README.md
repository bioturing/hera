<p align="center"><img src="https://user-images.githubusercontent.com/13636609/28245625-6356a0c2-6a35-11e7-8052-a936c1c928f2.png" width="50%"/></p>

Developed by BioTuring (www.bioturing.com), <i>hera</i> is a bioinformatics tool that helps analyze RNA-seq data. With a single command line, <i>hera</i> provides: 

- Base-to-base alignment BAM file
- Transcript abundance estimation
- Fusion gene detection with fused sequence assemblies 

Each process in <i>hera</i> was carefully organized and optimized in order to maximize the performance in term of time and accuracy. Hera quantification algorithm obtained the best ranking in a recent round of the SMC-RNA DREAM challenge: https://www.synapse.org/#!Synapse:syn2813589/wiki/423306 


# Example data
We designed a test using 20 datasets from Synapse Dream Challenge SMC-RNA, each of which contains 60 million read pairs. The test was done on a 32-core machine running Ubuntu 14.04. The result is shown in the table below:

<table width="100%">
   <tr>
      <td rowspan="11" width="400px">
         <img src="https://user-images.githubusercontent.com/13636609/28252091-a6d3126e-6ab6-11e7-90c4-2fee5f22716f.png" width="100%"/>
      </td>
      <td></td>
      <td align="center"> <b>Transcriptome</b> </td>
      <td align="center"> <b>Transcriptome + Genome</b> </td>
  </tr>
  <tr>
      <td colspan="3" align="center"> <b>Alignment</b> </td>
  </tr>
  <tr>
      <td align="center">Mapped read </td>
      <td align="center"> 93.3860% </td>
      <td align="center"> 93.3871%  </td>
  </tr>
  <tr>
      <td align="center">Memory </td>
      <td align="center"> 8GB </td>
      <td align="center"> 30GB </td>
  </tr> 
  <tr>
  <td colspan="3" align="center"> <b>Abundance estimation results</b> </td>
  </tr>
  <tr>
      <td align="center">Spearman</td>
      <td align="center">0.9033</td>
      <td align="center">0.9057</td>
  </tr>
  <tr>
      <td align="center">Pearson</td>
      <td align="center">0.9951</td>
      <td align="center">0.9951 </td>
  </tr>
  <tr>
  <td colspan="3" align="center"> <b>Gene fusion results</b> </td>
  </tr>
  <tr>
      <td align="center">True positive</td>
      <td align="center" colspan="2">0.6960</td>
  </tr>
  <tr>
      <td align="center">False negative</td>
      <td align="center" colspan="2">0.304</td>
  </tr>
  <tr>
      <td align="center">False positive</td>
      <td align="center" colspan="2">0.0595</td>
  </tr>
</table>

# Core algorithm

### Alignment
In <i>hera</i>, alignment starts with a hash-based approach that is applied on all the reads to anchor gene fragments. Then, Needleman–Wunsch algorithm is used to fill in the gaps between these anchored seeds. With this approach, the mapping time is reduced without hurting the precision. An additional conversion of transcriptome takes place to generate genome coordinates from the original transcriptome coordinates. This step provides a much better accuracy for splicing detection than mapping data onto a reference genome.

In another hand, <i>hera</i> is still able to perform the common genome mapping due to the incompletion of available transcriptomics. Any reads that cannot be mapped properly on the transcriptome will be remapped to the genome later. The procedure for this case is the same as the transcriptome mapping except the hash-based method is replaced with the Burrow-Wheeler Transform algorithm.

### Abundance estimation
Expectation–maximization algorithm is optimized with the SQUAREM procedure (Varadhan, R. & Roland, C. Scand. J. Stat. 35, 335–353 (2008)).
   
### Fusion detection
In order to detect fusions, <i>hera</i> keeps track of abnormally mapped reads. Based on their potential fusion site, these reads are divided into several groups, each of which is assembled into a super contig. These contigs will be mapped back onto the reference genome and thereby reveal their fusion gene pairs.

# Build requirements:
 
  * GNU GCC C Compiler
  * CMake (http://www.cmake.org/) version 3.1.0 or newer
  * liblzma-dev (Ubuntu) or xz-devel (Centos, Fedora, Red Hat) or xz (MacOS)
  * libbz2-dev (Ubuntu) or bzip2-devel.x86_64 (Centos, Fedora, Red Hat) or bzip2 (MacOS)
  * libz-dev (Ubuntu) or zlib-devel.x86_64 (Centos, Fedora, Red Hat) or zlib (MacOS)

# Install:

  ```shell
  1. git clone https://github.com/bioturing/hera.git
  2. cd hera/
  3. chmod +x build.sh
  4. ./build.sh
  ```

  For Linux users, <i> hera </i> can be easily installed by using <b>bioconda</b>:
  
  ```shell
  conda config --add channels bioconda
  conda install hera
  ```

# Usage:

### INDEX:
  ```
  ./hera/build/hera_build
          --fasta genome_sequence.fa (text file only)
          --gtf annotation_file.gtf
          --outdir path/to/output_directory
 [OPTIONAL]
          --full_index 0: None, 1: index full genome
          --grch38     0: No, 1: Yes
  ```
  
 By default, <i>hera</i> needs ~8GB for transcriptome indexing only. Full genome indexing needs ~30GB. You also can download indexed human genome file here: [GRCh37.75](https://drive.google.com/file/d/0B8iAyV1a-kmCRzZueVhGaDJuaDA/view), [GRCh38.82](https://drive.google.com/file/d/0ByiG7pl1z3EMT3F3UDBvR1JrR00/view)
 
 If you are running with GRCh38 human genome included ALT contigs, you should let hera know by defining --grch38 1.

### RUN:
  ```
  ./hera/build/hera quant -i path/to/index_directory [OPTIONAL] read1.fastq read2.fastq
  
  [OPTIONAL]:
    -o [output directory] (default: ./)
    -t [number of running threads] (default: 1)
    -z [level of bam file compression (1 - 9)] (default: -1)
    -b [Number of boostrap] (default: 1)
    -w [Output bam file 0: true, 1: false] (defaut: 0)
    -f [Genome fasta file]
   ```
  
  Eg: hera quant -i index/ -t 32 read1.fastq read2.fastq
  
  1. <b>Index directory</b>: Directory contain index file from previous index step
  
  2. <b>Genome fasta file</b>: If not defined, genome mapping will be ignore. Mapping on transcriptome needs ~8BG, but mapping with genome needs ~30GB.
  
  3. Output file include:
  - abundance.tsv  : Transcripts abundance estimation (tsv file)
  - abundance.h5   : Transcripts abundance estimation and boostrapping result (hdf5 file)
  - fusion.bedpe   : Fusion detection result (for paired-end data only)
  - transcript.bam : Alignment result

# Third-party

<i>hera</i> includes some third-patry software:
  * hdf5 [https://support.hdfgroup.org/HDF5/]
  * htslib [http://www.htslib.org/]
  * jemalloc [http://jemalloc.net/]
  * libdivsufsort [https://github.com/y-256/libdivsufsort]
  * zlib [https://zlib.net/]

# Contacts

Please report any issues directly to the github issue tracker. Also, you can send your feedback to info@bioturing.com

# Contributions
BioTuring Algorithm Team & 
Thao Truong, Khoa Nguyen, Tuan Tran, and Son Pham

# License

MIT license

Copyright (c) BioTuring Inc. 2017
All rights reserved. 
This Hera 1.0 version is freely accessible for both academic and industry users.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

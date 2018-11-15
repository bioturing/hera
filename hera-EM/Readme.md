This is Hera-EM code for RECOMB 2019.

It's important to note that Hera-EM module has not yet been incorporated into public Hera code yet.

While Hera is licensed under MIT, Hera-EM is not. To use the code, please send email to: info@bioturing.com

To compile code for linux, run sh ./build.sh

Transcriptome fasta for bowtie2 indexing can be download here: https://drive.google.com/file/d/10-weMz3BEShF45LgnAEDO88zfE_GHYck/view?usp=sharing

Running command: Hera-EM alignment.bam n_thread

alignment.bam is the output bam file from bowtie2 which is run as below:

bowtie2 -p n_thread -x index --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 --no-mixed --no-discordant -k 200 -I 1 -X 1000 -1 read1.fq -2 read2.fq | samtools view -bS - > alignment.bam

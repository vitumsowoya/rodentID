conda activate cytb_phylo
fastqc #to check quality of the file
bwa index cytb_db.fa
bwa mem -t 12 cytb_db.fa PE8L_S1_L001_R1_001.fastq.gz PE8L_S1_L001_R2_001.fastq.gz  | samtools view -b - | samtools sort -@ 12 -m 512M -o mapped.sorted.bam
samtools flagstat mapped.sorted.bam
samtools view -b -F 4 mapped.sorted.bam > mapped_only.bam
samtools flagstat mapped_only.bam
samtools sort mapped_only.bam -o mapped_only.sorted.bam
samtools index mapped_only.sorted.bam 
samtools idxstats mapped_only.sorted.bam > idxstats.txt
samtools fastq mapped_only.sorted.bam -1 mapped_R1.fastq.gz -2 mapped_R2.fastq.gz -0 /dev/null -s /dev/null -n
mkdir assembly
cd assembly
mv ../mapped_R1.fastq.gz . 
mv ../mapped_R2.fastq.gz . #move the mapped_R1.fastq.gz and mapped_R2.fastq.gz to the assembly folder
megahit -1 mapped_R1.fastq.gz -2 mapped_R2.fastq.gz -o megahit_out
cd megahit_out
cat final.contigs.fa
python3 ~/tools/quast-5.3.0/quast.py final.contigs.fa -o quast_output
mkdir blast 
cd blast
mv ../final.contigs.fa .
cp /home/vmsowoya/mbewa_study/Version_2/PE8L/cytb_3.fa .
makeblastdb -in cytb_3.fa  -dbtype nucl -out cytb_db
blastn   -query final.contigs.fa   -db cytb_db   -out cytb_vs_contigs.tsv   -outfmt 6   -evalue 1e-10   -num_threads 8
prodigal -i final.contigs.fa -a proteins.faa -d genes.fna -o prodigal.gff -p meta
mkdir phylogeny
cd phylogeny
cp ../final.contigs.fa .
seqkit grep -r -p "^k119_1\b" final.contigs.fa > PE8L_cytb_extract.fa #After this, open the new fa file and rename it by replacing contig name with rodent ID
cp /home/vmsowoya/mbewa_study/Version_2/PE8L/cytb_3.fa .
cat PE8L_cytb_extract.fa cytb_3.fa > cytb_4.fa
mafft --adjustdirection --auto cytb_4.fa > aligned.fa
trimal -in aligned.fa -out trimmed.fa
sed 's/^>_R_/>/' trimmed.fa > cleaned_4.fa
iqtree -s cleaned_4.fa -m MFP -B 1000 -nt AUTO


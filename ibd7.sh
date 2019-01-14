
f=$(sed -n "${SGE_TASK_ID}{p;q;}" 2017_gz.txt)
echo "THIS RUN IS--------------------------------- $f --------------"
DATA="/home3/redwards/IBD/HMP2/MGX/2017-12-14/"
mkdir temp_$f
newf=$(echo $f | sed -e s/.gz//g)
fastq=~/cpg/temp_$f/$newf
outtemp=~/cpg/temp_$f

## Unzip the file to a new temp directory
gunzip -c $DATA/$f > $fastq

## Use Perl to calculate kmers and substrings of kmers
#perl ~/cpg/IBD_Kmer_Analysis/GetNucFrequency_PerSeq_varK.pl $fastq 6 > $fastq.6mers
#~/bin/Python-3.6.0/python ~/cpg/IBD_Kmer_Analysis/kmerbias.py -f $fasta.6mers -k CGTACG  > $outtemp/bias.csv


## Use Python to count motifs
~/bin/Python-3.6.0/python IBD_Kmer_Analysis/count_motifs.py -i ~/cpg/temp_$f/*.fastq -o results/$f.motifs.csv

## Use Jellyfish to get raw kmers counts
/usr/local/jellyfish/bin/jellyfish count $fastq -C -m 1 -s 50M -o $outtemp/1mers.jf
/usr/local/jellyfish/bin/jellyfish count $fastq -C -m 2 -s 50M -o $outtemp/2mers.jf
/usr/local/jellyfish/bin/jellyfish count $fastq -C -m 3 -s 50M -o $outtemp/3mers.jf
/usr/local/jellyfish/bin/jellyfish count $fastq -C -m 4 -s 50M -o $outtemp/4mers.jf
/usr/local/jellyfish/bin/jellyfish dump $outtemp/1mers.jf > $outtemp/1mers.fasta                                                     
/usr/local/jellyfish/bin/jellyfish count $fastq -C -m 5 -s 50M -o $outtemp/5mers.jf
/usr/local/jellyfish/bin/jellyfish dump $outtemp/2mers.jf > $outtemp/2mers.fasta
/usr/local/jellyfish/bin/jellyfish count $fastq -C -m 6 -s 50M -o $outtemp/6mers.jf
/usr/local/jellyfish/bin/jellyfish dump $outtemp/3mers.jf > $outtemp/3mers.fasta                                                     
/usr/local/jellyfish/bin/jellyfish dump $outtemp/4mers.jf > $outtemp/4mers.fasta
/usr/local/jellyfish/bin/jellyfish dump $outtemp/5mers.jf > $outtemp/5mers.fasta
/usr/local/jellyfish/bin/jellyfish dump $outtemp/6mers.jf > $outtemp/6mers.fasta

# Count lines
wc -l $fastq > $outtemp/word_count.txt

$ clean up
rm $fastq

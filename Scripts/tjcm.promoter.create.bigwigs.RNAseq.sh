module load slurm
module load deeptools/3.3.1
module load WiggleTools/1.2.2

## mean of replicates, mean of KOs
## mean of reps
upperdir=/data/bio/tan/Mammalian2019Promoter #$1
# dir=CM2019RNAseqHuman #$2
# chromsize=/data/bio/tan/genome_fasta/genome.chrom.sizes

dir=CM2019RNAseqMouse #$2
chromsize=/data/bio/tan/genome_fasta/genome.chrom.sizes

cd $upperdir/$dir

##DEseq2 dexseq promoter counts 
tmp=deseq2.dexseq.promoter.txt

##DEseq2 dexseq exons counts 
tmp=deseq2.dexseq.counts.txt

##DEseq2 feature counts
tmp=deseq2.feature.counts.txt

##edger dexseq promoter counts 
tmp=edger.dexseq.promoter.txt

##edger dexseq exons counts 
tmp=edger.dexseq.counts.txt

##edger feature counts
tmp=edger.feature.counts.txt

##The rest steps (human use dexseq.counts, feature.counts; mouse use dexseq.counts, dexseq.promoter)
name=`echo $tmp |sed 's/.txt//g'`
mkdir -p bamCoverage/$name
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '2,50p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; sf=`echo $line |cut -d ":" -f2`; SlurmEasy -t 20 -d 264953 "bamCoverage --scaleFactor $sf -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '51,100p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; sf=`echo $line |cut -d ":" -f2`; SlurmEasy -t 20 -d 265003 "bamCoverage --scaleFactor $sf -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '101,150p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; sf=`echo $line |cut -d ":" -f2`; SlurmEasy -t 20 -d 265053 "bamCoverage --scaleFactor $sf -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '151,200p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; sf=`echo $line |cut -d ":" -f2`; SlurmEasy -t 20 -d 265103 "bamCoverage --scaleFactor $sf -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '201,250p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; sf=`echo $line |cut -d ":" -f2`; SlurmEasy -t 20 -d 265153 "bamCoverage --scaleFactor $sf -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '251,300p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; sf=`echo $line |cut -d ":" -f2`; SlurmEasy -t 20 -d 265203 "bamCoverage --scaleFactor $sf -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '301,318p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; sf=`echo $line |cut -d ":" -f2`; SlurmEasy -t 20 -d 265253 "bamCoverage --scaleFactor $sf -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
mkdir -p Bigwig_merged/$name
for samplename in `ls bamCoverage/$name/*bw |cat |cut -d "/" -f3 |cut -d_ -f1,2 |uniq`;
do SlurmEasy -t 1 "wiggletools mean bamCoverage/$name/${samplename}*bw | wigToBigWig stdin $chromsize Bigwig_merged/$name/${samplename}.bw";
done

##cpm normalization (for both human and mouse)
tmp=edger.feature.counts.txt
name=cpm.normalization
mkdir -p bamCoverage/$name
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '2,50p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; SlurmEasy -t 20 "bamCoverage --normalizeUsing CPM -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '51,100p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; SlurmEasy -t 20 -d 264349 "bamCoverage --normalizeUsing CPM -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '101,150p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; SlurmEasy -t 20 -d 264399 "bamCoverage --normalizeUsing CPM -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '151,200p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; SlurmEasy -t 20 -d 264449 "bamCoverage --normalizeUsing CPM -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '201,250p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; SlurmEasy -t 20 -d 264499 "bamCoverage --normalizeUsing CPM -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '251,300p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; SlurmEasy -t 20 -d 264549 "bamCoverage --normalizeUsing CPM -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '301,318p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; SlurmEasy -t 20 -d 264599 "bamCoverage --normalizeUsing CPM -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
mkdir -p Bigwig_merged/$name
for samplename in `ls bamCoverage/$name/*bw |cat |cut -d "/" -f3 |cut -d_ -f1,2 |uniq`;
do SlurmEasy -t 1 "wiggletools mean bamCoverage/$name/${samplename}*bw | wigToBigWig stdin $chromsize Bigwig_merged/$name/${samplename}.bw";
done


##Bigwig normalization using DESeq2 scale factor
# #SlurmEasy -t 100 "multiBamSummary bins --bamfiles STAR/*bam --scalingFactors bamCoverage/normal_normalization.txt -p 100 -o bamCoverage/normal_normalization.npz"
# #for samplename in `cut -f4 ProcessedData/PromoterObjects/human.sample.sheet.txt |sed '1d' |uniq`; do condition=`grep $samplename ProcessedData/PromoterObjects/human.sample.sheet.txt |cut -f1 |sed 's/$/.DESeq2.bw/g' |sed 's/^/bamCoverage\//g' |xargs`; echo "$samplename:$condition" >> bamCoverage/merge.txt; done && sed 's/ /;/g' bamCoverage/merge.txt > bamCoverage/merge.txt1 && mv bamCoverage/merge.txt1 bamCoverage/merge.txt
# tmp=multibamsummary_normalization.txt
# name=`echo $tmp |sed 's/.txt//g'`
# for line in `sed 's/\t/:/g' bamCoverage/$tmp| sed -n '11,50p'`; do input=`echo $line |cut -d ":" -f3 |sed 's/.bam//g'`; output=`echo $line |cut -d ":" -f1 |sed 's/.bam//g'`; sf=`echo $line |cut -d ":" -f2`; SlurmEasy -t 20 "bamCoverage --scaleFactor $sf -p 20 -b STAR/${input}.bam -o bamCoverage/$name/${output}.bw"; done
# ## merge replicates 
# mkdir Bigwig_merged/$name
# for i in `sed -n '1,30p' bamCoverage/merge.txt |sed 's/ /;/g'`;
# do condition=`echo $i |cut -d ":" -f1`; samplename=`echo $i |cut -d ":" -f2 |sed 's/;/ /g'`;
# SlurmEasy -t 1 "wiggletools mean ${samplename} | wigToBigWig stdin $chromsize Bigwig_merged/${condition}.bw";
# done



## pygenome Tracks
cd $upperdir/$dir
module load pyGenomeTracks/2.0 && module load bedtools2/2.27.0 && mkdir PygenomeSnapshot
bedtools sort -i /data/bio/tan/Ensembl/genes.bed > PygenomeSnapshot/genes.sorted.bed
bed=$upperdir/$dir/PygenomeSnapshot/genes.sorted.bed #/data/bio/tan/Ensembl/genes.bed #
rna_bw=$upperdir/$dir/Bigwig_merged/
make_tracks_file --trackFiles $rna_bw/*p25year* $bed -o $upperdir/$dir/PygenomeSnapshot/pygenome.tracks.ini
cd $upperdir/$dir/PygenomeSnapshot 
for i in `echo chr12:47657959-47708810`; do pyGenomeTracks --tracks pygenome.tracks.ini --region $i --width 50 -o RPAP3.gene.png; done

bed3=$upperdir/$dir/PygenomeSnapshot/brain.specific.bed
SlurmEasy -t 50 "pyGenomeTracks --tracks pygenome.tracks.ini --BED $bed3 --width 25 --outFileName brain25year"

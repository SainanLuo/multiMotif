iupac2meme -dna TGTTA >motif.me
fimo --verbosity 1 --thresh 1.0E-0 --oc fimo motif.me GCA_000814675.1.promoter.fa
cut -f 3,4,5,6,7,8,9,10 fimo/fimo.tsv|awk -F"\t" '{if($4~/+/ && $8~/TGTTA/){print $0}}' >fimo_out

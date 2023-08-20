echo -e "iupac2meme -dna TGTAAATTTACA >motif0.me" >>inpuc2meme.sh
for i in {1..40};do
	n_string=$(printf 'N%.0s' $(seq 1 $i))
	echo -e "iupac2meme -dna TGTAAA"$n_string"TTTACA >motif"$i".me" >>inpuc2meme.sh
done

sh inpuc2meme.sh

for i in {0..40};do
	fimo --verbosity 1 --thresh 1.0E-4 --oc fimo_"$i" motif"$i".me GCA_000008765.1.promoter.fa
done

cut -f 3,4,75,6,7,8,9,10 fimo_0/fimo.tsv |grep "TGTAAATTTACA" >>fimo_out

for i in {1..40};do
        n_string=$(printf 'N%.0s' $(seq 1 $i))
        cut -f 3,4,5,6,7,8,9,10 fimo_"$i"/fimo.tsv |awk -F"\t" '$8 ~ /^TGTAAA.*TTTACA$/ {print $0}' >>fimo_out
done

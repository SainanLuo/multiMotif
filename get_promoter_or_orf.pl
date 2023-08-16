#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($fna_file, $gff_file, $output_file, $upstream,$downstream, $promoter, $orf);

GetOptions(
	"fna=s"      => \$fna_file,      # 输入FNA文件
	"gff=s"      => \$gff_file,      # 输入GFF文件
	"up:i"        => \$upstream,       # 窗口长度（可选，默认值为400）
	"down:i"        => \$downstream,       # 窗口长度（可选，默认值为0）
	"o=s"        => \$output_file,   # 输出文件
	"promoter"   => \$promoter,      # 选择 promoter 参数
	"orf"        => \$orf,           # 选择 ORF 参数
	"h"          => \&help            # 帮助信息
);
# 打印帮助信息并退出
sub help {
	print "Usage: perl script.pl [options]\n";
	print "Options:\n";
	print "  -f, --fna      <file>    Input FNA file\n";
	print "  -g, --gff      <file>    Input GFF file\n";
	print "  -up, --upstream   <value>   gene start location upstream length (optional, default is 400)\n";
	print "  -down, --dowmstream   <value>   gene start location downstream length (optional, default is 0)\n";
	print "  -o, --output   <file>    Output file\n";
	print "  --promoter               Extract promoters\n";
	print "  --orf                    Extract ORFs\n";
	print "  -h                       Display this help and exit\n";
	exit;
}

# 检查必需的参数是否提供
unless ($fna_file && $gff_file && $output_file) {
	die "Usage: perl script.pl -f <fna_file> -g <gff_file> -o <output_file> [--promoter | --orf] [-l <windows>]\n";
}

open(Seq, $fna_file) or die "Cannot open $fna_file: $!";
open(IN, $gff_file) or die "Cannot open $gff_file: $!";
open(OUT, ">$output_file") or die "Cannot create $output_file: $!";
my %hash;
my $genome;

while (<Seq>) {
	chomp;
	if (/^>/) {
		my @array = split(" ", $_);
		$array[0] =~ s/>(.*)/$1/;
		$genome = $1;
		$hash{$genome} = '';
	} else {
		$hash{$genome} .= $_;
	}
}

$upstream ||= 400;  # 如果未提供窗口长度参数，则默认值为400
$downstream ||= 0;  # 如果未提供窗口长度参数，则默认值为0
$promoter = 1 unless ($promoter || $orf);  # 默认情况下选择 promoter

while (<IN>) {
	chomp;
	next if /^#/;
	my @temp = split /\t/, $_;
	my $name = $temp[0];
	my $left = $temp[3];
	my $right = $temp[4];
	my @arr = split(";", $temp[8]);
	my ($CdsID, $parent,$gene, $product) = ("None","None", "None", "None");
    
	foreach my $field (@arr) {
		if ($field =~ /^ID=(.*)/) {
			$CdsID = $1;
		}
		if ($field =~ /^Parent=(.*)/) {
                        $parent = $1;
                }
		if ($field =~ /^gene=(.*)/) {
			$gene = $1;
		} 
		if ($field =~ /^product=(.*)/) {
			$product = $1;
		}
	}
    
	if (exists $hash{$name}) {
		my $seqences="";
		if ($promoter) {
			if(exists $hash{$name}){
				if($temp[2] =~/CDS/ && $temp[6] eq "-"){
					my $start=$right - $downstream;
					my $windows=$upstream + $downstream;
					$seqences=substr($hash{$name},$start,$windows);
					$seqences=reverse($seqences);
					$seqences=~tr/ATCG/TAGC/;
					print OUT ">$name|$CdsID|$parent|$gene|$left|$right|$temp[6]|$product\n";
                                	print OUT "$seqences\n";
				}
				elsif($temp[2]=~/CDS/ && $temp[6] eq "+" ){
					my $start=$left-$upstream;
					if($start < 0){
						$start=0;
					}
					my $windows=$upstream + $downstream;
					$seqences=substr($hash{$name},$start,$windows);
					print OUT ">$name|$CdsID|$parent|$gene|$left|$right|$temp[6]|$product\n";
					print OUT "$seqences\n";
				}
				else{next;}
		
			}
		}
		elsif ($orf) {
			if(exists $hash{$name}){
				if($temp[2]=~/CDS/ && $temp[6] eq "-"){
					$seqences=substr($hash{$name},$left,$right-$left);
					$seqences=reverse($seqences);
					$seqences=~tr/ATCG/TAGC/;
					print OUT ">$name|$CdsID|$parent|$gene|$left|$right|$temp[6]|$product\n";
					print OUT "$seqences\n";
				}
				elsif($temp[2]=~/CDS/ && $temp[6] eq "+" ){
					$seqences=substr($hash{$name},$left,$right-$left);
					print OUT ">$name|$CdsID|$parent|$gene|$left|$right|$temp[6]|$product\n";
					print OUT "$seqences\n";

				}
				else{next;}
			}
      
		}
	}
}

close(IN);
close(Seq);
close(OUT);

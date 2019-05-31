#!/usr/bin/perl

# By I-Hsuan Lin
# Usage: perl run_trimming.pl folder_name raw 341F 805R
# fastq.gz are placed in folder_name/raw, all files are processed and trimmed files are placed in folder_name/trimmed

use strict;
use warnings;
use File::Basename;

my $path = shift @ARGV;
my $raw = shift @ARGV;
my $trimmed = "$path/trimmed";
my $tmp = "$path/tmp";

my $primerF = shift @ARGV;
my $primerR = shift @ARGV;
my $primerFrev;
my $primerRrev;

if($primerF eq "341F") {
	$primerF = "^CCTACGGRAGGCAGCAG";
	$primerFrev = "CTGCTGCCTYCCGTAGG";
} else {
	die "Invalid forward primer: $primerF\n";
}

if($primerR eq "805R") {
	$primerR = "^GACTACHVGGGTATCTAATCC";
	$primerRrev = "GGATTAGATACCCBDGTAGTC";
} else {
	die "Invalid reverse primer id: $primerR\n";
}

foreach my $f (glob("$path/$raw/*_R1_001.fastq.gz")) {

	my $file = basename($f);
	my $dir = dirname($f);
	my $oid = "";
	my $nid = "";
	my $lane = "";

	if($file !~ /^(\S+)_(S\d+_L\d+)_R1_001.fastq.gz/) {
		die "Invalid file pattern: $file\n";
	} else {
		$oid = $1;
		$lane = $2;
	}

	if(! -d "$trimmed") {
		print "mkdir $trimmed\n";
		mkdir "$trimmed";
	}

	if(! -d "$tmp") {
		print "mkdir $tmp\n";
		mkdir "$tmp";
	}

	$nid = ($oid =~ /^(\d)$/) ? "0$1" : $oid;

	my $fq1 = "$dir/$file";
	my $fq2 = "$dir/$oid\_$lane\_R2_001.fastq.gz";
	my $log = "$trimmed/$nid.log";
	my $opt = "";
	my $cmd = "";

	print "Running cutadapt on Sample $nid\n";

	$opt = "--match-read-wildcards --no-indels --error-rate 0.15 --overlap 8";	# Maximum error rate 0.15; Minimum overlap 8
	$cmd = "cutadapt $opt --discard-untrimmed -g $primerF -o $tmp/$nid.01.fastq.gz -p $tmp/$nid.02.fastq.gz $fq1 $fq2 > $log";
	print "$cmd\n";
	system($cmd);
	$cmd = "cutadapt $opt --discard-untrimmed -g $primerR -o $tmp/$nid.12.fastq.gz -p $tmp/$nid.11.fastq.gz $tmp/$nid.02.fastq.gz $tmp/$nid.01.fastq.gz >> $log";
	print("$cmd\n");
	system($cmd);

	$cmd = "cutadapt $opt -g $primerF -o $tmp/$nid.01.fastq.gz -p $tmp/$nid.02.fastq.gz $fq1 $fq2 > $log";
	print "$cmd\n";
	system($cmd);
	$cmd = "cutadapt $opt -g $primerR -o $tmp/$nid.12.fastq.gz -p $tmp/$nid.11.fastq.gz $tmp/$nid.02.fastq.gz $tmp/$nid.01.fastq.gz >> $log";
	print"$cmd\n";
	system($cmd);

	$opt = "--match-read-wildcards --no-indels --error-rate 0.15 --overlap 10 --minimum-length=200";
	$cmd = "cutadapt $opt -a $primerRrev -o $tmp/$nid.21.fastq.gz -p $tmp/$nid.22.fastq.gz $tmp/$nid.11.fastq.gz $tmp/$nid.12.fastq.gz >> $log";
	print "$cmd\n";
	system($cmd);
	$cmd = "cutadapt $opt -a $primerFrev -o $trimmed/$nid.2.fastq.gz -p $trimmed/$nid.1.fastq.gz $tmp/$nid.22.fastq.gz $tmp/$nid.21.fastq.gz >> $log";
	print"$cmd\n";
	system($cmd);
}

system("rm -Rf $tmp");

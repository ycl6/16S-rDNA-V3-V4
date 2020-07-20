#!/usr/bin/perl

# By I-Hsuan Lin
# Usage: perl run_trimming.pl project_folder fastq_folder forward_primer_sequence reverse_primer_sequence
#
# Example 1 - If "run_trimming.pl" is placed one-level up outside project_folder
# Example 1 directory structure:
# /home/user/my_ngs/run_trimming.pl			# Outside project_folder (default)
# /home/user/my_ngs/PRJEB27564				# A project
# /home/user/my_ngs/PRJEB27564/raw/*_1.fastq.gz		# raw fastq files, read 1
# /home/user/my_ngs/PRJEB27564/raw/*_2.fastq.gz		# raw fastq files, read 2
# /home/user/my_ngs/PRJEB27564/trimmed/*_1.fastq.gz	# cutadapt trimmed files, read 1
# /home/user/my_ngs/PRJEB27564/trimmed/*_2.fastq.gz	# cutadapt trimmed files, read 2
# Usage: perl run_trimming.pl PRJEB27564 raw CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC
#
# Example 2 - If "run_trimming.pl" is placed within the project_folder
# Example 2 directory structure:
# /home/user/my_ngs/PRJEB27564
# /home/user/my_ngs/PRJEB27564/run_trimming.pl		# Within project_folder
# /home/user/my_ngs/PRJEB27564/raw/*_1.fastq.gz
# /home/user/my_ngs/PRJEB27564/raw/*_2.fastq.gz
# /home/user/my_ngs/PRJEB27564/trimmed/*_1.fastq.gz
# /home/user/my_ngs/PRJEB27564/trimmed/*_2.fastq.gz
# Usage: perl run_trimming.pl . raw CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC
#

use strict;
use warnings;
use File::Basename;

my $path = shift @ARGV;
my $fastq = shift @ARGV;

$fastq = "$path/$fastq";
my $trimmed = "$path/trimmed";
my $tmp = "$path/tmp";

my $primerF = shift @ARGV;
my $primerR = shift @ARGV;
my $primerFrev = reverse $primerF;
$primerFrev =~ tr/acgtrymkbdhvACGTRYMKBDHV/tgcayrkmvhdbTGCAYRKMVHDB/;
my $primerRrev = reverse $primerR;
$primerRrev =~ tr/acgtrymkbdhvACGTRYMKBDHV/tgcayrkmvhdbTGCAYRKMVHDB/;

print "P1-FWD: $primerF\tP2-FWD: $primerR\n";
print "P1-REV: $primerFrev\tP2-REV: $primerRrev\n";

# Check Python version
my $pyj = "";
my $pyv = qx(python --version 2>&1);

if($pyv =~ /Python 3/) {
	print "Python 3\n";
	$pyj = "-j 0";	# allow auto-detection of CPU cores
} elsif($pyv !~ /Python 2/) {
	die "Command 'python' not found\n";
}

foreach my $f (glob("$fastq/*1.f*q.gz")) {
	my $file = basename($f);
	my $dir = dirname($f);
	my $id = "";
	my $fq1 = "$dir/$file";
	my $fq2 = "";

	# Accepts "_R1_001.fastq.gz" or "_1.fastq.gz"
	if($file =~ /^(\S+)_R1_001\.fastq\.gz/) {	# _R1_001.fastq.gz
		$id = $1;
		$fq2 = "$dir/$id\_R2_001.fastq.gz";
	} elsif($file =~ /^(\S+)_1\.fastq\.gz/) {	# _1.fastq.gz
		$id = $1;
		$fq2 = "$dir/$id\_2.fastq.gz";
	} else {
		die "Invalid file pattern: $file\n";
	}
	
	my $trim1 = "$trimmed/$id\_1.fastq.gz";
	my $trim2 = "$trimmed/$id\_2.fastq.gz";
	my $log = "$trimmed/$id.log";
	my $opt = "";
	my $cmd = "";

	if(! -d "$trimmed") {
		print "mkdir $trimmed\n";
		mkdir "$trimmed";
	}

	if(! -d "$tmp") {
		print "mkdir $tmp\n";
		mkdir "$tmp";
	}
	print "Running cutadapt on Sample \"$id\"\n";

	$opt = "-g $primerF -a $primerRrev -G $primerR -A $primerFrev -n 2 --match-read-wildcards --length 300 --minimum-length 150 --overlap 10 $pyj";
	$cmd = "cutadapt $opt -o $trim1 -p $trim2 $fq1 $fq2 > $log";
	print "$cmd\n";
	system($cmd);
}

system("rm -Rf $tmp");

#!/usr/bin/perl

# By I-Hsuan Lin
# Usage: perl lefse.pl lefse_table.res graphlan_outtree.txt expr.out

use strict;
use warnings;
use Data::Dumper;

if(!@ARGV) { die "Input is empty";}

my $res = shift @ARGV;  # lefse_table.res
my $txt = shift @ARGV;  # graphlan_outtree.txt
my $out = shift @ARGV;  # expr.out

my %resData;
my $i = 0;
open(RES, $res) or die "Cannot open $res";
while(<RES>) {
        chomp;
        my @data = split /\t/, $_;
        if($data[2] eq "") {
                next;
        }

        my $fulltaxon = $data[0];
        my $logmaxpct = $data[1];
        my $condition = $data[2];
        my $lda = $data[3];
        my $pvalue = $data[4];

        my @taxonomy = split /\./, $fulltaxon;
        my $rank = scalar(@taxonomy);
        my $last = pop @taxonomy;

        if($rank == 2) {
                $rank = "P";    # Phylum
        } elsif($rank == 3) {
                $rank = "C";    # Class
        } elsif($rank == 4) {
                $rank = "O";    # Order
        } elsif($rank == 5) {
                $rank = "F";    # Family
        } elsif($rank == 6) {
                $rank = "G";    # Genus
        }

	if(!defined($resData{$last}{'fulltaxon'})) {
		$i = 0;
	} else {
		$i++;
		$last = "$last.$i";
	}

        $resData{$last}{'fulltaxon'} = $fulltaxon;
       	$resData{$last}{'logmaxpct'} = $logmaxpct;
        $resData{$last}{'condition'} = $condition;
       	$resData{$last}{'lda'} = $lda;
        $resData{$last}{'pvalue'} = $pvalue;
       	$resData{$last}{'rank'} = $rank;
}
close RES;

my %txtData;
my $j = 0;
open(TXT, $txt) or die "Cannot open $txt";
while(<TXT>) {
        chomp;
        # <property applies_to="clade" datatype="xsd:string" id_ref="annotation" ref="A:1">A:Coriobacteriia</property>
        if($_ =~ /id_ref="annotation".+>(\w+):(.+)<\/property>/) {
                my $order = $1;
                my $name = $2;
                $name =~ s/ /_/g;
		if(!defined($txtData{$name})) {
			$j = 0;
		} else {
			$j++;
			$name = "$name.$j";
		}
                $txtData{$name} = $order;
        }
}
close TXT;

# Add alphabate id to %resData
my $resCount = 0;
my $txtCount = 0;
if(scalar(keys(%resData)) != scalar(keys(%txtData))) {
	# Fix potential Actinobacteria (Phylum/Class problem)
	$resCount = grep (/Actinobacteria/, keys(%resData));
	$txtCount = grep (/Actinobacteria/, keys(%txtData));

	if(!($resCount == 1 and $txtCount == 2)) {
		print "resData: ", scalar(keys(%resData)), "\ttxtData: ", scalar(keys(%txtData)), "\n";
        	die "Error: %resData != %txtData\n"
	}
}

foreach my $key (sort keys %txtData) {
	my $value = $txtData{$key};
	if($resCount == 1 and $txtCount == 2 and $key eq "Actinobacteria.1") {
		$key = "Actinobacteria";
	}

	if(defined($resData{$key}{'order'})) {
		$resData{$key}{'order'} .= ", $value";
	} else {
		$resData{$key}{'order'} = $value;
	}
}

# Print results
if(! defined $out or $out eq "") { die "Output filename not given\n" };
open(OUT, ">$out") or die "Cannot open $out\n";

print OUT "fulltaxon\ttaxon\torder\trank\tcondition\tlda\tpvalue\tlogmaxpct\n";
foreach my $taxon (keys %resData) {
	my $new_taxon = $taxon;
	$new_taxon =~ s/\.\d$//;
        print OUT "$resData{$taxon}{'fulltaxon'}\t$new_taxon\t$resData{$taxon}{'order'}\t$resData{$taxon}{'rank'}\t$resData{$taxon}{'condition'}\t";
        print OUT "$resData{$taxon}{'lda'}\t$resData{$taxon}{'pvalue'}\t$resData{$taxon}{'logmaxpct'}\n";
}
close OUT;

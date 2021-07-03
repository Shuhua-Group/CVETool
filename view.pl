#!/usr/bin/perl
use List::Util qw/max min/;
use strict;

my $loc=$ARGV[0];
my $infile=$ARGV[1];
my @SV;
my @SNP;
my @exp;
my %expn;
my @out;

open(INFILE,$infile);
while(<INFILE>){
	push @SV,(split /\t/)[0];
	push @SNP,(split /\t/)[1];
	push @exp,(split /\t/)[4];
}
chomp @exp;
close INFILE;

open(OUT,">>view.tmp0");

my $svmax = max @SV;
my $snpmax = max @SNP;

for (my $i = 0; $i <= $svmax; $i++) {
	for (my $j = $snpmax; $j >= 0; $j--) {
		my $var="SV${i}_SNV${j}";
		$expn{"$var"}=$var;
		for(my $index=1;$index<@SV;$index++){
			if(@SV[$index] == $i and @SNP[$index] == $j){
				$expn{"$var"}=$expn{"$var"}."\t".@exp[$index];
			}
		}
	}
}
my $maxcount=0;
for (my $i = 0; $i <= $svmax; $i++) {
	for (my $j = $snpmax; $j >= 0; $j--) {
		my $var="SV${i}_SNV${j}";
		my @tmp=split("\t",$expn{"$var"});
		my $sampcount=@tmp-1;
		if ($sampcount > 0) {
			print OUT $var."_".$sampcount."\t";
		}
		if ($sampcount > $maxcount) {
			$maxcount = $sampcount;
		}
	}
}
print OUT "\n";

for (my $i = 0; $i <= $svmax; $i++) {
	for (my $j = $snpmax; $j >= 0; $j--) {
		my $var="SV${i}_SNV${j}";
		my @tmp=split("\t",$expn{"$var"});
		my $sampcount=@tmp-1;
		if ($sampcount > 0) {
			for (my $l = 0; $l < $sampcount; $l++) {
				@out[$l]=@out[$l]."\t".@tmp[$l+1];
			}
			for (my $l = $sampcount; $l < $maxcount; $l++) {
				@out[$l]=@out[$l]."\tNA";
			}
		}
	}
}
for(my $i=0;$i<@out;$i++){
	print OUT @out[$i]."\n";
}
close (OUT);
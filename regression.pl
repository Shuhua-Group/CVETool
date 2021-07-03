#!/usr/bin/perl
use strict;
my $loc=$ARGV[0];
my $snpfile=$ARGV[1];
my $svfile=$ARGV[2];
my $expfile=$ARGV[3];
my $covfile=$ARGV[4];
my $rscript=$ARGV[5];
my $Rsoft=$ARGV[6];
my %gtsnp;
my %gtsv;
my %exp;
my %gender;
my %age;

open(INFILE,"${loc}/${snpfile}");
my @snpfile=<INFILE>;
chomp @snpfile;
close INFILE;
my @snphead=split("\t",@snpfile[0]);
my @snpbody=split("\t",@snpfile[1]);
for (my $i = 2; $i < @snphead; $i++) {
	$gtsnp{"@snphead[$i]"}=@snpbody[$i];
}

open(INFILE,"${loc}/${svfile}");
my @svfile=<INFILE>;
chomp @svfile;
close INFILE;
my @svhead=split("\t",@svfile[0]);
my @svbody=split("\t",@svfile[1]);
for (my $i = 3; $i < @svhead; $i++) {
	$gtsv{"@svhead[$i]"}=@svbody[$i];
}

open(INFILE,"${expfile}");
my @expfile=<INFILE>;
chomp @expfile;
close INFILE;
my @exphead=split("\t",@expfile[0]);
my @expbody=split("\t",@expfile[1]);
for (my $i = 1; $i < @exphead; $i++) {
	$exp{"@exphead[$i]"}=@expbody[$i];
}

open(INFILE,"${covfile}");
my @covfile=<INFILE>;
chomp @covfile;
close INFILE;
my @covhead=split("\t",@covfile[0]);
my @genderbody=split("\t",@covfile[1]);
my @agebody=split("\t",@covfile[2]);
for (my $i = 1; $i < @covhead; $i++) {
	$gender{"@covhead[$i]"}=@genderbody[$i];
	$age{"@covhead[$i]"}=@agebody[$i];
}

open(OUT,">>${loc}/SV_SNP.paired.txt");
print OUT "#CHR\tSNPpos\tSVstart\tSVend\tgene\tregresssion.p\tSV.coef\tSV.p\tSNP.coef\tSNP.p\n";

open(TMP,">>${loc}/SV_SNP.paired.tmp");
print TMP "SV\tSNP\tgender\tage\texp\n";

for(my $i=2;$i<@snphead;$i++){
	my $cursamp=@snphead[$i];
	if(!$gender{"$cursamp"} and $gender{"$cursamp"} ne "0" or !$age{"$cursamp"} and $age{"$cursamp"} ne "0" or !$exp{"$cursamp"} and $exp{"$cursamp"} ne "0"){
		next;
	}else{
		print TMP $gtsv{"$cursamp"}."\t".$gtsnp{"$cursamp"}."\t".$gender{"$cursamp"}."\t".$age{"$cursamp"}."\t".$exp{"$cursamp"}."\n";
	}
}
close (TMP);

my $tmpout=`Rscript ${rscript} ${loc}/SV_SNP.paired.tmp`;
my @out=split(" ",$tmpout);
chomp @out;
print OUT @snpbody[0]."\t".@snpbody[1]."\t".@svbody[1]."\t".@svbody[2]."\t".@expbody[0]."\t".@out[0]."\t".@out[1]."\t".@out[2]."\t".@out[3]."\t".@out[4]."\n";

close(OUT);

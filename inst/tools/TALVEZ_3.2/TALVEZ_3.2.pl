#!/usr/bin/perl -w

################################################################
# TALVEZ 3.2  							#
# by Alvaro L Perez-Quintero                       		#
#   < alperezqui at gmail dot com >                		#
#   Created 6 - VI - 2011 					#
#   Last modification 10 - VII - 2012      			#
#								#
#This script uses java code from Matzieu and Hatzigeorgiu 2010  #
#available at:							#
#http://www.diana.pcbi.upenn.edu/Tools/PlantTFBS/  		#
################################################################

#this version doesn't calculate minimum thresholds for each TAL but rather keeps the same user-defined scoring threshold
#this version uses RVD-DNA matrices supplied by the user as separate files
#this version searches both strands for EBEs

use warnings;
use strict;
use Data::Dumper;
#use Fatal qw(open close);
use Getopt::Std;

######################################  INPUTS

#optional parameters

my $minimumscore;
my $pseudocounts;
my $minimumlength;
my $numbertargets;
my $matrix_file1;
my $matrix_file2;
my $execution_dir;

my %options=();
getopts("ht:a:l:v:e:z:p:", \%options);

if (defined $options{t}) {$numbertargets = $options{t}} else {$numbertargets = 0};
if (defined $options{a}) {$minimumscore = $options{a}} else {$minimumscore = 6};
if (defined $options{l}) {$minimumlength = $options{l}} else {$minimumlength = 0};
if (defined $options{v}) {$pseudocounts = $options{v}} else {$pseudocounts = 0.0001};
if (defined $options{e}) {$matrix_file1 = $options{e}} else {$matrix_file1 = 'mat1'};
if (defined $options{z}) {$matrix_file2 = $options{z}} else {$matrix_file2 = 'mat2'};
if (defined $options{p}) {$execution_dir = $options{p}} else {$execution_dir = 'tmp'};

my $IN0;
my $IN1;

my $help = "\nCorrect usage: perl TALVEZ_2.5pl (-t -a -l -v -e -z -p) TALfile Promoters(Fasta) (output)\n
TALfile: tab separated >Talid\tRVDs(XX-XX-XX)
Promoters: Fasta file, sequences in one line

Parameters:

-t : Number of Targets to output for each TAL, ranked according to score (default = all)
-a : minimum score Allowed to report binding sites (default = 6)
-l : length RVDs before position correction, initial repetition to apply a scaled-down matrix to acount for mismatches in the N' terminal region (default = 0)
-v : pseudocount number (binding probability for any base regardless of the RVD) incorporated to Positional weight matrix (default = 0.0001)
-e : RVD-DNA specificities to use for RVDs before position correction, format (tab -separated) RVD A C G T, (default = mat1)
-z : RVD-DNA specificities to use for RVDs after position correction, format (tab -separated) RVD A C G T, (default = mat2)
-p : project name to be used as working directory (default = tmp)
\n\n";



#####open files

# Data file with TAL sequences, each sequence is in the format 'NN-NI-HD..', each sequence has a FASTA-like header, header and sequence are tab separated, with a .txt extension. IMPORTANT: there should be no spaces in the header names

if ($ARGV[0]) {open $IN0, "< $ARGV[0]" or die "Can not open $ARGV[0]!!\n";} 
	 else {print "$help";
		exit;} 

# FASTA file with promoter sequences. IMPORTANT: fasta sequences should comprise just one line

if ($ARGV[1]) {open  $IN1, "< $ARGV[1]"or die "Can not open $ARGV[1]!!\n";}
	 else {print $help;
		exit;}

if (!-d $execution_dir)
{
	mkdir("$execution_dir");
}

my $nuname = $ARGV[0];
$nuname =~ s/\.txt//g;
open my $OUT0, ">$execution_dir/$nuname"."_pwms.txt"; #files for java scripts

print "-t number of Targets to print # $numbertargets\n";
print "-a minumum Allowed score $minimumscore\n";
print "-l Length (RVDs) used for position correction  $minimumlength\n";
print "-v pseudocounts $pseudocounts\n";
print "-e Matrix1 $matrix_file1\n";
print "-z Matrix2 $matrix_file2\n";



####################  RVD-DNA CODE TO USE for PMWs ############

#Reads a file with the RVD-DNA specificities in the format (tab -separated) "RVD A C G T" and asigns them to a hash %RVDs

open my $MAT1, "< $matrix_file1" or die "Can not open $matrix_file1\n";
my %RVDs;

while (my $line = <$MAT1>){
        chomp $line;
	my @tmparray = split ("\t", $line);
	$RVDs{$tmparray[0]} = join ("\t", $tmparray[1],$tmparray[2],$tmparray[3],$tmparray[4]);
			}
close $MAT1;


######################### RVD-DNA CODE Scaled-down TO USE for position correction

open my $MAT2, "< $matrix_file2" or die "Can not open $matrix_file2\n";

my %PCRVDs;

while (my $line = <$MAT2>){
        chomp $line;
	my @tmparray = split ("\t", $line);
	$PCRVDs{$tmparray[0]} = join ("\t", $tmparray[1],$tmparray[2],$tmparray[3],$tmparray[4]);
			}
close $MAT2;

#variables for PMWs
my $amin;

############################################ CREATE PWMS #################################################################

while(my $line = <$IN0>){
	chomp $line;
	my @line1 = split (/\t/, $line); # Get TAL sequences
	print $OUT0 $line1[0]."\n";
	my @TALRVDs = split (/\-/, $line1[1]); #get RVDs
	print $OUT0 "=\t$RVDs{OO}\n"; #First repetition T or C
		my $nulength;
		if ($#TALRVDs > $minimumlength and $minimumlength > 0) {$nulength = $minimumlength -1;} else {$nulength = $#TALRVDs};
		for my $j (0 .. $nulength) {
		$amin = $TALRVDs[$j];
		if (defined $RVDs{$amin}) {print $OUT0 "=\t$RVDs{$amin}\n"} else {print $OUT0 "=\t$RVDs{XX}\n"};
			}
		for my $j ($nulength+1 .. $#TALRVDs) {   #use a scaled-down version of the RVD code matrix toprint PMWs
		$amin = $TALRVDs[$j];
		if (defined $PCRVDs{$amin}) {print $OUT0 "=\t$PCRVDs{$amin}\n"} else {print $OUT0 "=\t$RVDs{XX}\n"};
			}
	print $OUT0 "\n";
}

close $IN0;
close $OUT0;

#################### Check fasta format and create array from FASTA file to find and extract locations in parser

my @array_lines_fasta;
my @sequence;
my $seq;
my @header;

open my $Revcom, ">$execution_dir/RevComp_$ARGV[1]"; #Reverse_complemented file

while (my $line = <$IN1>){
        chomp $line;
        push @array_lines_fasta, $line;
			}

close $IN1;

my %fullheader;

for (my $i = 0; $i <= $#array_lines_fasta; $i += 2)
	{
        push @sequence, $array_lines_fasta[$i+1];

	my @tmpheader = split (" ", $array_lines_fasta[$i]);

	my $bla = $tmpheader[0];
	$bla =~ s/^>//g;
	$fullheader{$bla} = "$array_lines_fasta[$i]";
	$fullheader{"$bla"."_revcom"} = "$array_lines_fasta[$i]";

	push @header , $tmpheader[0];
	

	print $Revcom "$tmpheader[0]\n";
	print $Revcom "$array_lines_fasta[$i+1]\n";

	$seq = $array_lines_fasta[$i+1];
	$seq =~ tr/atcgATCG/tagcTAGC/;
 	$seq = reverse($seq);
	
	print $Revcom "$tmpheader[0]_revcom\n";	
	print $Revcom "$seq\n";

	push @header , "$tmpheader[0]_revcom";
	
	push @sequence, $seq;
	if ($tmpheader[0] =~ /^>/) {next;} else {die "Fasta format is incorrect, plase re-format to have sequences in just one line";}

#revcom
	}

close $Revcom;

my %index_header = map {$header[$_] => $_ } 0..$#header; #create index of headers


#################### Find tal lengths and create array from TAL file for parser

open $IN0, "< $ARGV[0]";
my @sequencetal;
my @headertal;
my @TALlengths;

while (my $line = <$IN0>){
        chomp $line;
	my @tmparray = split ("\t", $line);
        push @headertal, $tmparray[0];
	push @sequencetal, $tmparray[1];
			}

close $IN0;

my %index_headertal = map {$headertal[$_] => $_ } 0..$#headertal; #create index of tal headers

foreach (@sequencetal) {my $seq = $_;
			$seq =~ s/-//g;
			my $length = (length $seq)/2;
			push @TALlengths, $length;
			}

####################################### Run PlantTFBS java scripts ############################################

my $PMWs = "$execution_dir/$nuname"."_pwms.txt";
my $Back = "$execution_dir/$ARGV[1]"."background";
my $Thresh = "$execution_dir/$nuname"."ThresholdsTALs";


system ("java -Xms96m -Xmx2048m simplescancode.Background $execution_dir/RevComp_$ARGV[1] $Back");


########################### added by alexis ################################
system("sed 's/\,/\./g' $Back >$Back.2");
rename("$Back.2",$Back);
############################################################################


#create threshold file
open my $OUT2, "> $Thresh";                                                    
foreach  (@headertal){
			$_ =~ s/^>//g;
			print $OUT2 "$_\t$minimumscore\n";
			}
close $OUT2;

###################SCREEN FOR BINDING SITES

my $OUTPUTname;
if ($ARGV[2]) {$OUTPUTname = $ARGV[2];} else {$OUTPUTname = "output"};


system ("java -Xms96m -Xms2048m simplescancode.Scan $pseudocounts $PMWs $Thresh $Back $execution_dir/RevComp_$ARGV[1] $execution_dir/$OUTPUTname");

########################### added by alexis ################################
system("sed 's/\,/\./g' $execution_dir/$OUTPUTname.scores >$execution_dir/$OUTPUTname.scores.2");
rename("$execution_dir/$OUTPUTname.scores.2","$execution_dir/$OUTPUTname.scores");
############################################################################

################################################### PARSE OUTPUT   #############################################

open my $IN3, "< $execution_dir/$OUTPUTname.scores"; #file scores
open my $IN4, "< $execution_dir/$OUTPUTname.locations"; #file locations
open my $OUT3, "> $execution_dir/$OUTPUTname"."arranged.tmp";
open my $OUT4, "> $OUTPUTname"."_complete"; #output

my $talname;
my @TALNAME;	
my @arrangedcolumns;
my @arrangedlocations;


while (my $line = <$IN3>)
{
        chomp $line; #removes \n
	if ($line =~ />/)  {$talname = $line; #asigns lines starting with ">"
				@TALNAME = split ("\t", $talname); 
				}		
	else 		{
			my $line1 = "$TALNAME[0]$TALNAME[1]\t$line"; #appends the tal name to lines not starting with '>'
       				my @arraycolumns = split ("\t", $line1);
			for (my $i = 2; $i <= $#arraycolumns; $i++) #asigns number to the i counter from 2 to column numbers
				{
				push  @arrangedcolumns, "$arraycolumns[0]\t$arraycolumns[1]\t$arraycolumns[$i]";
				} #prints each column of a line after the first column
			}
}
while (my $line = <$IN4>)
{
        chomp $line; #removes \n
	if ($line =~ />/) {next;} 
		else    {
			my @arraycolumns = split ("\t", $line);
			for (my $i = 1; $i <= $#arraycolumns; $i++)
			{push @arrangedlocations, "$arraycolumns[$i]\n";}
			}
}

my @Unsorted;

for (my $i = 0; $i <= $#arrangedcolumns; $i++){
#print $OUT3 "$arrangedcolumns[$i]\t$arrangedlocations[$i]";
push @Unsorted, "$arrangedcolumns[$i]\t$arrangedlocations[$i]";
						}

my @sorted=sort {(split(/\t/,$a))[0] cmp (split(/\t/,$b))[0] || (split(/\t/,$b))[2]<=>(split(/\t/,$a))[2]} @Unsorted;

my $rank ;

for (my $i = 0; $i <= $#sorted; $i++)
			{
		chomp $sorted[$i];
		my @a = split(/\t/,"$sorted[$i]");
		my @b = split(/\t/,"$sorted[$i-1]");
			if ("$a[0]" eq "$b[0]") {$rank++}
			else {$rank = 1}
		print $OUT3 "$sorted[$i]\t$rank\n";
			}


close $OUT3;
close $IN3;
close $IN4;

#search for the target sequence in the fasta and for the tal length in tal array

my $chromosome; #fasta header of the sequence that includes the target sequence
my @arrayTalvez; #Tal name \t target header (chromosome) \t score \t location
my $TALID; #Tal header
my $score;
my $seq_array_number; #variable to identify index number of header in FASTA array
my $tal_array_number;  #variable to identify index number of header in TAL array
my $seqtoscan; #sequence containing target sequence
my $tallength;
my $targetseq;
my $inipos; #initial psition of the TAL boz, counting from the end of the sequence to the beggining
my $TALSEQ;

my $geneinfo;

###############Get Annotations for genes from annotation files corresponding to the scanned files


my $Annotfile =$ARGV[1];
my @tmpname = split ("_", $Annotfile);
$Annotfile = "$tmpname[0]";
$Annotfile = "$Annotfile.Annotations";

my %Annotations;

if (-e $Annotfile) {open my $IN6, "$Annotfile";
while (my $line = <$IN6>){
        chomp $line;
        my @tmparray = split ("\t", $line);
        my @tmparray2;
        for my $i (1 .. $#tmparray) {push @tmparray2, $tmparray[$i]};
        $Annotations{$tmparray[0]} = join ("\t", @tmparray2);
                        }
close $IN6;}
else {print "Gene Annotations not available\n";}

##################Final parsed file


open my $IN5, "$execution_dir/$OUTPUTname"."arranged.tmp";
	
print $OUT4 "TAL_ID\tTAL_SEQ\tSEQ_ID\tSCORE\tEBEstrand\tTALBS_distance_from_end\tTALBS_start\tTALBS_end\tTALBS_sequence\tRANK\tChromosome\tStrand\tGene_start\tGene_end\tAnnotation\n";

my $chromosome2; #full ID including spaces
my $rinipos;
my $rinipos1;
my $rfinpos;

while (my $line = <$IN5>)
{	
	chomp $line;
	@arrayTalvez = split ("\t", $line);
	$TALID = "$arrayTalvez[0]";
	$chromosome = "$arrayTalvez[1]";
	$chromosome2 = $fullheader{"$arrayTalvez[1]"};
	$score = "$arrayTalvez[2]";
	$inipos ="$arrayTalvez[3]";
	$rank = "$arrayTalvez[4]";
	$tal_array_number = $index_headertal{"$TALID"};
	$tallength = $TALlengths["$tal_array_number"];

if ($chromosome =~ "_revcom") {


	$TALSEQ = $sequencetal["$tal_array_number"];
		$seq_array_number = $index_header{">$chromosome"}; #get position of chromosome in array of FASTA
		$seqtoscan = $sequence[$seq_array_number]; #find sequence in fasta
		my $inipos1 = (length $seqtoscan) + $inipos;
		$rinipos = ((length $seqtoscan) + $inipos)*-1; #distance in the positive strand
		$rfinpos = (length $seqtoscan) + $rinipos;
		my $finpos = $inipos1 + $tallength + 1;
		$rinipos1 = $rfinpos - $tallength -1;
		$targetseq = substr ($seqtoscan, $inipos, $tallength + 1); #extract target sequence
		if (defined $Annotations{$chromosome2}) {$geneinfo = $Annotations{$chromosome2};} else {$geneinfo = "-\t-\t-\t-\t-";}
	if ($numbertargets == 0) {
	print $OUT4 "$TALID\t$TALSEQ\t$chromosome2\t$score\t-strand\t$rinipos\t$rinipos1\t$rfinpos\t$targetseq\t$rank\t$geneinfo\n";
				}
	elsif ($rank <= $numbertargets) {
	print $OUT4 "$TALID\t$TALSEQ\t$chromosome2\t$score\t-strand\t$rinipos\t$rinipos1\t$rfinpos\t$targetseq\t$rank\t$geneinfo\n";
				}
	else {next;}


				} else

				{
	$TALSEQ = $sequencetal["$tal_array_number"];
		$seq_array_number = $index_header{">$chromosome"}; #get position of chromosome in array of FASTA
		$seqtoscan = $sequence[$seq_array_number]; #find sequence in fasta
		my $inipos1 = (length $seqtoscan) + $inipos;
		my $finpos = $inipos1 + $tallength + 1;
		$targetseq = substr ($seqtoscan, $inipos, $tallength + 1); #extract target sequence
		if (defined $Annotations{$chromosome2}) {$geneinfo = $Annotations{$chromosome2};} else {$geneinfo = "-\t-\t-\t-\t-";}
	if ($numbertargets == 0) {
	print $OUT4 "$TALID\t$TALSEQ\t$chromosome2\t$score\t+strand\t$inipos\t$inipos1\t$finpos\t$targetseq\t$rank\t$geneinfo\n";
				}
	elsif ($rank <= $numbertargets) {
	print $OUT4 "$TALID\t$TALSEQ\t$chromosome2\t$score\t+strand\t$inipos\t$inipos1\t$finpos\t$targetseq\t$rank\t$geneinfo\n";
				}
	else {next;}


				}
}
close $OUT4;
close $IN5;
exit;

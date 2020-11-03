#!/usr/bin/perl -w
use warnings;
use strict;
use Statistics;
use Getopt::Std;
use Statistics::R;
use List::Util qw( min max );
use List::MoreUtils qw(uniq);
use Bio::Perl;

##################################################################
# QueTAL - FuncTAL v1.1						 #
# by Alvaro L Perez-Quintero                       		 #
#   < alperezqui at gmail dot com >                		 #
#   Created:17/4/2015						 #
#  					     			 #
#								 #
#This script uses subroutines from the program compareMotifs.pl  #
#from the Homer suite available at http://biowhat.ucsd.edu/homer/#
##################################################################


#Choose scoring method for similarity matrix, correlation recommended (Pearson)

#$scoreMethod = 'absDifference';
#$scoreMethod =  'freqError';
my $scoreMethod = 'correlation';
my $norevopp = 1; #option to analyze similarities in the opsite binding site, 1 = only same sense

#Help info

sub printCMD() {
print STDERR "\n

############################################################################
			QueTAL - FuncTAL v1.1: 

This script takes a series of RVD sequences from TAL effectors 
and calculates distances based on DNA binding specificities

This script uses subroutines from the program compareMotifs.pl 
from the Homer suite available at http://biowhat.ucsd.edu/homer/

Requires: R, Perl module Statistics::R, R library 'ape'
#############################################################################

Usage: FuncTAL.pl [options] TALfile Outname(default= output) AdditionalInfo(optional)

TALfile: Fasta file containing nucleotide or protein sequences for TAL effectors 
OR tab separated  file containing one TAL sequence perline in format: >Talid\tHG-NG-NN...\tAdditional_info_column(optional)
AdditionalInfo: Tab separated file with fasta Ids in the first column and information 
to be used to color the output tree in the second column (e.g. specie or strain)


-n :[f, d, u, p or c] Neighbor_joining tree format: organze the output tree in a fan (f), phylogra, (p), dendrogram (d), unrooted (u)(default) or cladogram (c) format
-c :[T or F]: compare the input TALs against a set of available TAL sequences from the ncbi database
-l :[integer]: Limit, if c=T, number of target TAls form database to output (default = 5)

Outputs:
Outname.matrix: Distance matrix based on TAL repeats
Outname.tre: nj tree in newick format
Outname.png: nj tree in png format
Outname.hits (if -c = T): Similar TALs found in the database according to specificities 


";
exit;
}

if (@ARGV < 1) {
	printCMD();
}


#######Get arguments

my %options=();
getopts("hf:n:c:l:", \%options);

my $Treetype = '';
my $compare = '';
my $RVDmatrix = '';
my $lim = 5; #number of similar TALs to output if compare option is activated

if (defined $options{f}) {$RVDmatrix = $options{f}} else {$RVDmatrix = 'Info/2014mat18'};
if (defined $options{n}) {$Treetype  = $options{n}} else {$Treetype = "u"};
if (defined $options{c} && ("$ARGV[0]" ne 'Info/PublicTALs.txt')) {$compare = $options{c}} else {$compare = "F"};
if (defined $options{l}) {$lim = $options{l}} else {$lim = 5};

my $cmpMatrixFile;
my $Outdir = "Outputs";
if ( -d "$Outdir") {} else { system ("mkdir $Outdir")};
if ($ARGV[1]) {$cmpMatrixFile = "$Outdir/$ARGV[1]"} else {$cmpMatrixFile = "$Outdir/Output"};

####################Detect input file and format as needed
my @INPUT;

open my $IN0, "<$ARGV[0]" or die "Can not open $ARGV[0], file doesn't exist\n";


my @RAW;
my $tabfileflag = '';

while ( my $line = <$IN0> ) 	{
	chomp $line;
	push @RAW, $line;
		}


###Detect input format

my $countaa = 0;
my $countnuc = 0;
my $countfas = 0; 
for (my $i = 0; $i <= 10; $i++) {
 if ((defined $RAW[$i])  and ($RAW[$i] =~ '>')) {$countfas++;}
	elsif ((defined $RAW[$i]) and ("$RAW[$i]" =~ /[HDLRP]+/)) {$countaa++}
	elsif ((defined $RAW[$i]) and ("$RAW[$i]" =~ /[ATCG]+/i)) {$countnuc++}

}

#check if file already contains RVDs

if ($countfas == 0) {die "Input file doesn't contain fasta headers \n";}
elsif ($countnuc == 0 && $countaa == 0) { 
	for (my $i = 0; $i <= 10; $i++) {
	if  (defined $RAW[$i]) {
	my @tmparr = split ("\t", "$RAW[$i]");
	if ((defined $tmparr[1]) and ("$tmparr[1]" =~ /[NH-]+/)) {$tabfileflag = "T"}	
	
}#if defined
}#for
}#elsif

if ($tabfileflag eq "T") {@INPUT = @RAW; print STDERR "\tInput file is tabular file containing RVDs\n"}


####If file is fasta


else {

#convert to aminoacid if necessary

my %AA;
my $AAref; my $IDsref;

if ($countaa <= $countnuc) {
print STDERR "Input file is a nucleotide fasta, converting to aminoacid and getting RVDs\n";
($AAref, $IDsref) = read_fasta(\@RAW, "translate");
}#if nucleotide

else 

{
print STDERR "Input file is an aminoacid fasta, getting RVDs\n";
($AAref, $IDsref) = read_fasta(\@RAW);
}

%AA = %$AAref;
my @rIDs = @$IDsref;

open my $Inital, "<Info/Initial_TAL_strings.txt";
my @Ini_repeat;

while ( my $line = <$Inital> ) {
	chomp $line;
	push @Ini_repeat, $line;}

####################FORMAT REPEATS

for ( my $i = 0 ; $i <= $#rIDs ; $i++ ) {
my $id = "$rIDs[$i]";
my $seq = $AA{$rIDs[$i]};

my @position_array; #array that will store all the postitions where a string for a TAL repetition is found
#look for each possible TAL repetition initiation
foreach (@Ini_repeat) {
		my $ini_motif = $_;
		my $indexresult = 0;
		my $offset = 0;
  #use index function in a loop to find all occurences of the TAL initial string
		while ( $indexresult != -1 ) {
			$indexresult = index( $seq, $ini_motif, $offset );			
			$offset = $indexresult + 1;			
			if ( $indexresult > -1 ) { push @position_array, $indexresult }
		}
	}
  #create array with uniq positions where the the initial TAL string was found
	my @POS = uniq( sort { $a <=> $b } (@position_array) );

  #get repeats lengths
	my @Repeat_len;

	for ( my $i = 1 ; $i <= $#POS ; $i++ ) {
		my $length = $POS[$i] - $POS[ $i - 1 ];
		push @Repeat_len, $length;
		
	}
  #find most common repetition length
		my(%count);
		foreach my $value (@Repeat_len) {
    		$count{$value}++;
		}
		my $common_len = (sort {$count{$b} <=> $count{$a}} @Repeat_len)[0];
		my $min_len = min @Repeat_len;
		my $max_len = max @Repeat_len;

   #exclude sequences for which no TAL repetition was found 
		if (not defined $min_len) 
	{print "no TAL repetitions were found for sequence $id, it was excluded from further analyses\n"}
		else {

	my @RVD_array;
	my $rvd;
	for ( my $i = 0 ; $i <= $#Repeat_len ; $i++ ) {
		my $times = ($Repeat_len[$i] / $common_len);
		if ($times > 1.5) { #deal with XXXX repeats
		my $ini = "$POS[$i]";
		for  (my $k = 1 ; $k <= $times ; $k++ ) {
						$rvd = substr ($seq, $ini + 11, 2); 
						$ini = $ini + $common_len;
						push @RVD_array, $rvd;
							}}
		elsif ($Repeat_len[$i] == $common_len - 1) {$rvd = substr ($seq, $POS[$i] + 11, 1); push @RVD_array, $rvd."*"} 
		else {
		$rvd = substr ($seq, $POS[$i] + 11, 2);
		push @RVD_array, $rvd;}
						}
	#get last repeat		
		my $last_rvd = substr ($seq, $POS[$#POS] + 11, 2);
		push @RVD_array, $last_rvd;

	#print repeats
	
	#push @INPUT, "$id\t"; foreach (@RVD_array) {push @INPUT, "$_-"}; push @INPUT, "\n";

	#print "$id\t"; for ( my $l = 0 ; $l <= $#RVD_array - 1; $l++ ){print "$RVD_array[$i]-"}; print "$RVD_array[$#RVD_array ]\n";
	
	$INPUT[$i] = "$id\t"; 
	for ( my $l = 0 ; $l <= $#RVD_array - 1; $l++ ){$INPUT[$i] = join ("", $INPUT[$i], "$RVD_array[$l]-")}; 
	$INPUT[$i] = join ("", $INPUT[$i], "$RVD_array[$#RVD_array ]\n");
	print "$INPUT[$i]";
				}



}


}#else fasta
 
########### RVD-DNA CODE TO USE for PMWs ##############

#Reads a file with the RVD-DNA specificities in the format (tab -separated) "RVD A C G T" and asigns them to a hash %RVDs

open my $MAT1, "< $RVDmatrix" or die "Can not open $RVDmatrix, RVD-DNA specificities file doesn't exist\n";
my %RVDs;

while (my $line = <$MAT1>){
        chomp $line;
	my @tmparray = split ("\t", $line);
	$RVDs{$tmparray[0]} = join ("\t", $tmparray[1],$tmparray[2],$tmparray[3],$tmparray[4]);
			}
close $MAT1;

############## CREATE PWMS #####################################


my @PMW;
my $Addinfoflag = "F";
my @Addinfo;
my @OIDs;
my %hash; #to look for Ids for compare

print STDERR "\tCreating positional weight matrices...\n";

foreach my $line (@INPUT){
	chomp $line;
	my @line1 = split (/\t/, $line); # Get TAL sequences
	push (@PMW, $line1[0]."\n");
	$line1[0] =~ s/>//g;
	push (@OIDs, "$line1[0]");
	@hash{@OIDs}=();
	if ($compare eq "T") {push (@Addinfo, "0_UserInput"); $Addinfoflag = "T"}
	elsif ((defined $line1[2]) && ($compare eq "F")) {push (@Addinfo, "$line1[2]")}
	my @TALRVDs = split (/\-/, $line1[1]); #get RVDs
	if ((not defined $TALRVDs[1]) && (length $line1[1] > 2)) {die " Format error, make sure input file is in the format >Talid\tHG-NG-NN...)";}
		for my $j (0 .. $#TALRVDs) {
		my $amin = $TALRVDs[$j];
		if (defined $RVDs{$amin}) {push (@PMW, "$RVDs{$amin}\n")} else {push (@PMW, "$RVDs{XX}\n")};
			}
}
close $IN0;

if ($#Addinfo == $#INPUT) {$Addinfoflag = "T"};

my %pubRVDs;
my %pubSpecie;

if ($compare eq "T") {
open OUTcomp, ">$cmpMatrixFile.hits" or die "Couldn't open output file $cmpMatrixFile.hits\n";
open my $IN1, "<Info/PublicTALs.txt"; #file containing publicly availabe TALs to use for comparison
@hash{@OIDs}=(); #Aditional hash to identify the users original TALs
while(my $line = <$IN1>){
	chomp $line;
	my @line1 = split (/\t/, $line); # Get TAL sequences
	push (@PMW, $line1[0]."\n");
	$line1[0] =~ s/>//g;
	push (@Addinfo, "$line1[2]");
	$pubSpecie{$line1[0]} = "$line1[2]\t$line1[3]";
	$pubRVDs{$line1[0]} = $line1[1];
	my @TALRVDs = split (/\-/, $line1[1]); #get RVDs
		for my $j (0 .. $#TALRVDs) {
		my $amin = $TALRVDs[$j];
		if (defined $RVDs{$amin}) {push (@PMW, "$RVDs{$amin}\n")} else {push (@PMW, "$RVDs{XX}\n")};
			}
}
}

#Read additional Info file!!!!!

if ((defined $ARGV[2]) && ($compare eq "F") && ($Addinfoflag eq "F")) {
my $count = 0;
open my $ADD, "<$ARGV[2]" or die "Can not open $ARGV[0], file doesn't exist\n";
while ( my $line = <$ADD> ) {
	chomp $line;
	my $addID = (split (/\t/, $line))[0];
	$addID =~ s/>//g;
	if (exists $hash{$addID}) {
	$line =~ s/>//g;
	push @Addinfo, (split (/\t/, $line))[1];
	$count ++;
		}

}
if ($count == keys %hash) {$Addinfoflag = "T"} else {print STDERR "\tIDs in additional info file do no match input fasta file, additional info was not used\n"}
}


#####COMPARE MOTIFS

my $rawMotifs = readMotifFile(@PMW);

print STDERR "\tDetermining similar motifs...\n";
compMotifs($rawMotifs, $cmpMatrixFile);

#####BUILD TREE using R

print STDERR "\tBulding tree...\n";

my $R = Statistics::R->new(); #bridge to use with R

my $cmds1 = <<EOF;
library("ape")
exprs <- as.matrix(read.table("$cmpMatrixFile.matrix", header=TRUE, sep = "\t",row.names = 1, as.is=TRUE))
tr <- bionj(exprs)
#bt<-boot.phylo(tr, exprs, FUN = function(xx) bionj(exprs), B = 1000)
write.tree(tr, file="$cmpMatrixFile.tre")
EOF

my $run1 = $R->run($cmds1);

if ($Addinfoflag eq 'T') { #Use additional info to color the tree

$R->set('labels', \@Addinfo);

my $cmds2 = <<EOF;
flabels<-as.factor(labels)
colvector<-rainbow (nlevels(flabels))[as.integer(flabels)]
png("$cmpMatrixFile.png")
plot(tr,"$Treetype", tip.color=colvector, show.tip.label=TRUE, cex=0.6, use.edge.length = TRUE)
add.scale.bar()
#nodelabels(bt,frame = "none", cex = 0.5, col = "blue")
legend(x = 'bottomright',legend = as.character(levels(flabels)),col = rainbow (nlevels(flabels)), pch = par("pch"), bty = 'n', xjust = 1, cex = 0.7)
dev.off()
EOF
my $run2 = $R->run($cmds2);
}


else {

my $cmds2 = <<EOF;
png("$cmpMatrixFile.png")
plot(tr,"$Treetype", show.tip.label=TRUE, cex=1, use.edge.length = TRUE)
add.scale.bar()
#nodelabels(bt,frame = "none", cex = 0.5, col = "blue")
dev.off()
EOF
my $run2 = $R->run($cmds2);
}

#print "$out2";

############# Read motif File

sub readMotifFile {
	#my ($file) = @_;
	my @array = @_;
	my @motifs = ();

	#open IN, $file or die "Could not open file: \"$file\"\n";
	my $m='';
	my $count = 0;
	my $badFlag = 0;
	#while (<IN>) {
	foreach (@array) {
		$count++;
		chomp;
		s/\r//g;
		my @line= split /\t/;
		if ($line[0] =~ s/^>//) {
			if ($count > 1) {
			push(@motifs, $m);
			}
			$badFlag= 0;
			my @matrix = ();
			$m={matrix=>\@matrix};
			$m->{'cons'} = '';
			$m->{'name'} = $line[0];
			next;
		} else {
			if (@line != 4) {
				print STDERR "Wrong file format\n";
				exit;
			}
			push(@{$m->{'matrix'}}, \@line);
			$m->{'len'}++;
		}

	}
	if (scalar($m->{'matrix'}) > 0) {  #adds last motif?
				push(@motifs, $m);
					}
	#close IN;

	open OUTC, ">$cmpMatrixFile.cons" or die "Couldn't open output file $cmpMatrixFile.cons\n";

	for (my $i=0;$i<@motifs;$i++) {
		cleanUpMotif($motifs[$i]);
	print OUTC "$motifs[$i]->{'name'}\t$motifs[$i]->{'cons'}\n";
	}
	return \@motifs;

}


############ Create consesus sequence

sub cleanUpMotif {
	my ($motif) = @_;
	my $consensus = '';
	for(my $i=0;$i<$motif->{'len'};$i++) {
		my $sum = 0;
		my $code = 'A';
		my $best = 0;
		for (my $j=0;$j<4;$j++) {
			my $v = $motif->{'matrix'}->[$i][$j];
			$sum+= $v;
			if ($v > $best) {
				$best = $v;
				$code = 'A' if ($j==0);
				$code = 'C' if ($j==1);
				$code = 'G' if ($j==2);
				$code = 'T' if ($j==3);
			}
			if ($best < 0.4) {
				$code = 'N';
			}
		}
		for (my $j=0;$j<4;$j++) {
			$motif->{'matrix'}->[$i][$j] /= $sum;
		}
		$consensus .= $code;
	}
	$motif->{'cons'} = $consensus;
	if ($motif->{'name'} eq '') {
	$motif->{'name'} = $consensus;


	}


}


#######compMotifs

sub compMotifs {
	my ($motifs, $cmpMatrixFile) = @_;
	my @m = @$motifs;
	my $startNum = scalar(@m); #number of motifs
	my @newList = ();
	my $lastNUM = -1;
	my @cmpMatrix = ();
	if ($cmpMatrixFile ne '') {
		for (my $i=0;$i<@m;$i++) {
			my @a = ();
			for (my $j=0;$j<@m;$j++) {
				if ($i==$j) {
					push(@a, 1);
				} else {
					push(@a, 0);
				}
			}
			
			push(@cmpMatrix, \@a);
		}
	}

	for (my $i=0;$i<@m;$i++) {
		if ($i==$lastNUM) {
			print STDERR "BLAH = $i\n";
		}
		$lastNUM = $i;
		push(@newList, $m[$i]);
		for (my $j=$i+1;$j<@m;$j++) {
			my ($s,$o,$d,$L) = (0,0,0,0);
				($s,$o,$d,$L) = compareMotifs($m[$i], $m[$j]);
			if ($scoreMethod ne 'correlation') {
				$s = $s/$L;
			}
			if ($cmpMatrixFile ne '') {
				$cmpMatrix[$i][$j] = $s;
				$cmpMatrix[$j][$i] = $s;
			}

		}
	}

	if ($cmpMatrixFile ne '') {
		open OUT, ">$cmpMatrixFile.matrix" or die "Couldn't open matrix output file $cmpMatrixFile.matrix\n";
		print OUT "Matrix";
		for (my $i=0;$i<@m;$i++) {
			my $n = $m[$i]->{'name'};
			print OUT "\t$n";
		}
		print OUT "\n";
		for (my $i=0;$i<@m;$i++) {
			my $n = $m[$i]->{'name'};
			print OUT "$n";
			my %scores;
			my %cons;
			for (my $j=0;$j<@m;$j++) {
				#my $inv = 1 - $cmpMatrix[$i][$j];
				my $inv = 1 - $cmpMatrix[$i][$j];
					my $o = $m[$j]->{'name'};
					my $p = $m[$j]->{'cons'};
				print OUT "\t$inv";
				if ($compare eq "T" && exists $hash{$n} && not exists $hash{$o}) {$scores{$o} = $inv; $cons{$o}=$p;}
						}
			if ($compare eq "T" && exists $hash{$n}){ #find closest TALs in database
			print OUTcomp "The TALs with more similar binding sites to $n :".$m[$i]->{'cons'}." in our database are:\n";
			print OUTcomp "ID(genebank)\tDistance\tbindingsite\tRVDs\tSpecie\tStrain\n";
			my @SimilarID; my @SimilarScore;
			foreach my $name (sort { $scores{$a} <=> $scores{$b} } keys %scores) {push @SimilarID, "$name"; push @SimilarScore, "$scores{$name}";}
			for (my $j=0;$j<=$lim;$j++)
						{
	print OUTcomp "$SimilarID[$j]\t$SimilarScore[$j]\t$cons{$SimilarID[$j]}\t$pubRVDs{$SimilarID[$j]}\t$pubSpecie{$SimilarID[$j]}\n";
						}
								}
			print OUT "\n";
		}
		close OUT;
	}
		
}

#Compare motifs

sub compareMotifs  {
	my ($m1, $m2) = @_;

	my $bestOffset = 0;
	my $bestDirection = 0;
	my $bestScore = -1e10;

	my $rv2 = revoppMotif($m2);

	my @default = (0,0,0,0);
	my $len1 = $m1->{'len'};
	my $len2 = $m2->{'len'};
	my $offset2 = $len2-1;
	my $offset1 = 0;
	my $bestLength = 0;
	for (my $offset2 = $len1-1;$offset2>-1*$len2;$offset2--) {
		my $officialOffset = -1*$offset2+$offset1;
		#print STDERR "$officialOffset\n";
		my $max1 = $len1 - $offset1;
		my $max2 = $len2 - $offset2;
		my $curLen = $max1;
		if ($max2 < $curLen) {
			$curLen = $max2;
		}
		my @mm1 = ();
		for (my $i=$offset2;$i<0;$i++) {
			push(@mm1, \@default);
		}
		for (my $i=0;$i<$len1;$i++) {
			push(@mm1, $m1->{'matrix'}->[$i]);
		}
		for (my $i=0;$i<$offset2+$len2-$len1;$i++) {
			push(@mm1, \@default);
		}
		my @mm2 = ();
		my @mm2r = ();
		for (my $i=0;$i<$offset2;$i++) {
			push(@mm2, \@default);
			push(@mm2r, \@default);
		}
		for (my $i=0;$i<$len2;$i++) {
			push(@mm2, $m2->{'matrix'}->[$i]);
			push(@mm2r, $rv2->{'matrix'}->[$i]);
		}
		while (scalar(@mm2) < scalar(@mm1)) {
			push(@mm2, \@default);
			push(@mm2r, \@default);
		}
		my $curLength = @mm1;


		my $n2 = scalar(@mm2);
		my $n1 = scalar(@mm1);;

		my $score = scoreComparison(\@mm1, \@mm2);
		
		#print STDERR "$offset2 ($len1)\t$offset2 ($len2)\t$officialOffset\t$score\n";
		
			if ($score > $bestScore) {
				$bestScore = $score;
				$bestOffset = $officialOffset;
				$bestDirection = 0;
				$bestLength = $curLength;
			}
		
	}
	#print STDERR "$bestScore\t$bestOffset\t$bestDirection\t$bestLength\n";
	$bestOffset *= -1;
	return ($bestScore, $bestOffset, $bestDirection,$bestLength);
}

sub printMatrix {
	my ($m) = @_;
	print STDERR "Matrix:\n";
	foreach(@$m) {
		foreach(@$_) {
			print STDERR "\t$_";
		}
		print STDERR "\n";
	} 
}

sub scoreComparison {
	my ($m1, $m2) = @_;
	my $len = @$m1;
	my $score = 0;
	if ($scoreMethod eq 'absDifference') {
		for (my $i=0;$i<$len;$i++) {
			for (my $j=0;$j<4;$j++){ 
				$score += $m1->[$i][$j] * $m2->[$i][$j] - 0.0625;
			}
		}
	} elsif ($scoreMethod eq 'freqError') {
		for (my $i=0;$i<$len;$i++) {
			my $curScore = 0;
			my $expectedScore = 0;
			for (my $j=0;$j<4;$j++) { 
				my $diff = $m1->[$i][$j]-$m2->[$i][$j];
				$curScore -= $diff*$diff;
				for (my $k=0;$k<4;$k++) { 
					my $diff2 = $m1->[$i][$j]-$m2->[$i][$k];
					$expectedScore -= $diff2*$diff2/4;
				}
			}
			$score +=  -1*($curScore-$expectedScore)/$expectedScore;
		}
	} elsif ($scoreMethod eq 'correlation') {
		my @a1 = ();
		my @a2 = ();
		foreach(@$m1) {
			foreach(@$_) {
				push(@a1,$_);
			}
		}
		foreach(@$m2) {
			foreach(@$_) {
				push(@a2,$_);
			}
		}
		my $lp = 0;
		($score,$lp) = Statistics::correlation(\@a1,\@a2);
	}
	return $score;
}

sub revoppMotif {
	my ($m) = @_;
	my @a = ();
	my $rv = {matrix=>\@a,name=>$m->{'name'},cons=>$m->{'cons'},
					v=>$m->{'v'},p=>$m->{'p'},logp=>$m->{'logp'},str=>$m->{'str'},
					gapinfo=>$m->{'gapinfo'}, len=>$m->{'len'},mask=>$m->{'mask'},
					homer2=>$m->{'homer2'}};
	for (my $i=-1+$m->{'len'};$i>=0;$i--){ 
		my @b = ();
		for (my $j=3;$j>=0;$j--) {
			push(@b, $m->{'matrix'}->[$i][$j]);
		}
		push(@{$rv->{'matrix'}}, \@b);
	}
	revoppConsensus($rv);

	return $rv;
}
sub revoppConsensus {
	my ($motif) = @_;
	my $new = reverse($motif->{'cons'});
	$new =~ s/A/X/g;
	$new =~ s/T/A/g;
	$new =~ s/X/T/g;

	$new =~ s/C/X/g;
	$new =~ s/G/C/g;
	$new =~ s/X/G/g;

	$new =~ s/R/X/g;
	$new =~ s/Y/R/g;
	$new =~ s/X/Y/g;

	$new =~ s/M/X/g;
	$new =~ s/K/M/g;
	$new =~ s/X/K/g;

	$new =~ s/B/X/g;
	$new =~ s/V/B/g;
	$new =~ s/X/V/g;

	$new =~ s/D/X/g;
	$new =~ s/H/D/g;
	$new =~ s/X/H/g;
	$motif->{'cons'} = $new;
}

sub correlation {
	my ($x,$y) = @_;

	my $n = scalar(@$x);
	return 0 if ($n==0 || scalar(@$y) == 0);
	my $xysum = 0;
	my $xsum = 0;
	my $ysum = 0;
	my $x2sum = 0;
	my $y2sum = 0;
	for (my $i=0;$i<$n;$i++) {
		$xysum += $x->[$i]*$y->[$i];
		$xsum += $x->[$i];
		$ysum += $y->[$i];
		$x2sum += $x->[$i]*$x->[$i];
		$y2sum += $y->[$i]*$y->[$i];
	}
	my $r = "NA";
	my $numerator = $xysum - $xsum*$ysum/$n;
	my $denomerator = ($x2sum - $xsum*$xsum/$n)*($y2sum-$ysum*$ysum/$n);
	return ("NA","NA") if ($denomerator <= 0);
	$denomerator = sqrt($denomerator);
	if ($denomerator > 0) {
		$r = $numerator / $denomerator;
	}
	my $df = $n-2;
	my $den = ((1-$r)*(1+$r));
	if ($den < 1e-10) {
		$den = 1e-10;
	}
	my $t = $r*sqrt($df/($den));
	my $logp = 1;
	if (abs($t) > 1e-5) {
		$logp = logbetai(0.5*$df,0.5,$df/($df+$t*$t));
	}

	return ($r,$logp);

}




###reading fasta



sub read_fasta { #takes two arguments; one fasta array and a flag to translate or not to protein
my $inlines = shift;
my @lin= @{$inlines};
my @rIDs;
my $trans = shift;
my %fast;
for (my $i = 0; $i <= $#lin - 1; $i++)
			{
if ($lin[$i] =~ '>') {#process fasta
	push @rIDs, "$lin[$i]";
my $j = $i + 1;
my $string = '';
do
				{
		if ((defined $lin[$j]) && ($lin[$j] !~ '>')) {
			$string = join ("", $string, "$lin[$j]");
				$j++;
					}
else {$j= 0}
		}while($j != 0);
if (defined $trans) {my $prot_obj; $prot_obj = translate_as_string ($string); $fast{$lin[$i]} = "$prot_obj";}
else {$fast{$lin[$i]} = "$string"};
}
}
return (\%fast, \@rIDs);
}

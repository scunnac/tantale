#!/usr/bin/perl -w
use warnings;
use strict;
use Statistics::R;
use Getopt::Std;
use Statistics::Basic qw(:all);


##################################################################
# Evolve.pl							 #
# by Alvaro L Perez-Quintero                       		 #
#   < alperezqui at gmail dot com >                		 #
#   Created:17/4/2015						 #
#  					     			 #
##################################################################


sub printCMD() {
print STDERR "\n

Usage: Evolve Repeatfile(Info/ALLUniqaa)

############################################################################
			Evolve.pl

This scripts simulates TAL effector evolution by creating artificial TALs form
a collection of TAL repeats. It feed the descendent TALs into the programs DisTAL
and ClustalW. It reauires DisTAL ClustalW and all their dependencies


-n :[number] Number of trees to generate for each alpha and beta combination (default = 50)
-i :[number]: indel penalization for DisTAL/ARLEM (default = 10)
-d :[number]: duplication penalization for DisTAL/ARLEM (default = 10)
-l :[number]: length of initial TAL effectors in number of repeats (defaule = 15)

#############################################################################

";
exit;
}

if (@ARGV < 1) {
	printCMD();
}


open my $OUT3, ">Simulated_evolution_results";
#print $OUT3 "alpha\tbeta\tgamma\tdistmean\tdistmedian\tdiststdev\tprogram\n";

for (my $x = 0; $x <= 100; $x+=5) { #alues for alpha

my $probmut = $x / 1000;

for (my $y = 0; $y <= 100; $y+=5) {#values for beta

my $probdel = $y / 1000;

my @arrayclust;
my @arrayquetal;

my %options=();
getopts("hn:i:d:l:", \%options);

my $Iter = '';
my $indel = '';
my $dup = '';
my $len = '';

if (defined $options{n}) {$Iter = $options{n}} else {$Iter = 10}; #Number of Simulated TALs to get average for each alpha and beta combination
if (defined $options{i}) {$indel = $options{i}} else {$indel = 10};
if (defined $options{d}) {$dup = $options{d}} else {$dup = 10};
if (defined $options{l}) {$len = $options{l}} else {$len = 15};

for (my $s = 1; $s <= $Iter; $s++) {

my $TALlength = $len -1;

my $mutcicles = 100; #cycles of mutation per generation

my %Repeats;



my $maxwidth = ($TALlength/2); #allowed meximum length in repeats of an insterion, duplication or deletion in a cycle

open my $IN0, "<$ARGV[0]" or die "Can not open $ARGV[0], file doesn't exist\n";

while ( my $line = <$IN0> ) 	{
	chomp $line;
	my @tmp = split (/\t/, $line);
	$Repeats{$tmp[0]}=$tmp[1];
		}

my $count = keys %Repeats;

my @ancestor;

#Create ancestor sequences of length TALlength by selecting random numbers from a collection of uniq repeats

for (my $i = 1; $i <= $TALlength; $i++) {
	my $rep = int (rand ($count));
	push @ancestor, $rep;
	#print "$rep\t";
	}




#Mutate, GENERATE TWO DESCENDENTS

#Array of Array

my @AoA;

for my $i (0..1) {
	my @array = @ancestor;
	$AoA[$i] = [ @array ];
		  }



for (my $j = 0; $j <= $#AoA; $j++) { #for each descendent

	for (my $i = 0; $i <= $mutcicles; $i++){ #for each generation
		my $len = $#{$AoA[$j]};

		if ($len > 2) {

		my $randmut = rand (1);

		if ($randmut <= $probmut) 	{ #mutation
			my $randrepeat = int (rand ($len)); #choose a repeat to mutate
			my $rep = int (rand ($count)); #choose a new repat to mutate it to
			$AoA[$j][$randrepeat] = "$rep";
						}

		#my $randup = rand (1);

		#if ($randup <= $probdup) 	{ #duplication
		#	my $randrepeat = int (rand ($len)); #choose a repeat to duplicate
		#	my $randwith = int (rand ($maxwidth));
		#	for my $i (0 .. $randwith) { 
		#	splice @{$AoC[$j]}, $randrepeat, 0, "$AoC[$j][$randrepeat+$randwith]";
		#				}
		#				}



		my $randel = rand (1);

		if ($randel <= $probdel) 	{ #deletion
			my $randrepeat = int (rand ($len)); #choose a start point for deletion
			my $randwith = int (rand ($maxwidth));
			my $randinspoint = int (rand ($len)); #choose an insertion point
			my $randevent = rand (1);
			if ($randevent <= 0.5) {splice @{$AoA[$j]}, $randrepeat, $randwith}
			else {
			for my $i (0 .. $randwith) {splice @{$AoA[$j]}, $randinspoint, 0, "$AoA[$j][$randrepeat+$randwith]";}
			     }

						}

			}
	}
						
}



#Mutate, GENERATE FOUR DESCENDENTS of the original 2


my @AoB;

for my $i (0..1) {
	my @array = @{$AoA[0]};
	$AoB[$i] = [ @array ];
		  }
for my $i (2..3) {
	my @array = @{$AoA[1]};
	$AoB[$i] = [ @array ];
		  }




for (my $j = 0; $j <= $#AoB; $j++) { #for each descendent

	for (my $i = 0; $i <= $mutcicles; $i++){ #for each generation
		
		my $len = $#{$AoB[$j]};

		if ($len > 3) {

		my $randmut = rand (1);

		if ($randmut <= $probmut) 	{ #mutation
			my $randrepeat = int (rand ($len)); #choose a repeat to mutate
			my $rep = int (rand ($count)); #choose a new repat to mutate it to
			$AoB[$j][$randrepeat] = "$rep";
						}

		#my $randup = rand (1);

		#if ($randup <= $probdup) 	{ #duplication
		#	my $randrepeat = int (rand ($len)); #choose a repeat to duplicate
		#	my $randwith = int (rand ($maxwidth));
		#	for my $i (0 .. $randwith) { 
		#	splice @{$AoC[$j]}, $randrepeat, 0, "$AoC[$j][$randrepeat+$randwith]";
		#				}
		#				}




		my $randel = rand (1);

		if ($randel <= $probdel) 	{ #deletion
			my $randrepeat = int (rand ($len)); #choose a start point for deletion
			my $randwith = int (rand ($maxwidth));
			my $randinspoint = int (rand ($len)); #choose an insertion point
			my $randevent = rand (1);
			if ($randevent <= 0.5) {splice @{$AoB[$j]}, $randrepeat, $randwith}
			else {
			for my $i (0 .. $randwith) {splice @{$AoB[$j]}, $randinspoint, 0, "$AoB[$j][$randrepeat+$randwith]";}
			     }

						}
		}
	}
						
}


#Mutate, GENERATE EIGHT DESCENDENTS of the original 2


my @AoC;

for my $i (0..1) {my @array = @{$AoB[0]}; $AoC[$i] = [ @array ];}
for my $i (2..3) {my @array = @{$AoB[1]}; $AoC[$i] = [ @array ];}
for my $i (4..5) {my @array = @{$AoB[2]}; $AoC[$i] = [ @array ];}
for my $i (6..7) {my @array = @{$AoB[3]}; $AoC[$i] = [ @array ];}


for (my $j = 0; $j <= $#AoC; $j++) { #for each descendent

	for (my $i = 0; $i <= $mutcicles; $i++){ #for each generation
		
		my $len = $#{$AoC[$j]};

		if ($len > 2) {

		my $randmut = rand (1);

		if ($randmut <= $probmut) 	{ #mutation
			my $randrepeat = int (rand ($len)); #choose a repeat to mutate
			my $rep = int (rand ($count)); #choose a new repat to mutate it to
			$AoC[$j][$randrepeat] = "$rep";
						}

		#my $randup = rand (1);

		#if ($randup <= $probdup) 	{ #duplication
		#	my $randrepeat = int (rand ($len)); #choose a repeat to duplicate
		#	my $randwith = int (rand ($maxwidth));
		#	for my $i (0 .. $randwith) { 
		#	splice @{$AoC[$j]}, $randrepeat, 0, "$AoC[$j][$randrepeat+$randwith]";
		#				}
		#				}


		my $randel = rand (1);

		if ($randel <= $probdel) 	{ #deletion
			my $randrepeat = int (rand ($len)); #choose a start point for deletion
			my $randwith = int (rand ($maxwidth));
			my $randinspoint = int (rand ($len)); #choose an insertion point
			my $randevent = rand (1);
			if ($randevent <= 0.5) {splice @{$AoC[$j]}, $randrepeat, $randwith}
			else {
			for my $i (0 .. $randwith) {splice @{$AoC[$j]}, $randinspoint, 0, "$AoC[$j][$randrepeat+$randwith]";}
			     }

						}
						
}

}


open my $OUT0, ">Evolved";


##################################Get Sequences

my @Names = qw(A B C D E F G H);


for (my $j = 0; $j <= $#AoC; $j++) {
	my @array = @{$AoC[$j]};
	print $OUT0 ">$Names[$j]\n";
	foreach (@array) {print $OUT0 "$Repeats{$_}\n"}
}

for my $i ( 0 .. $#AoC ) {
     print STDERR "\t@{$AoC[$i]}\n";
 }


#########Run DisTAL

system ("perl DisTAL_v1.0.pl -m T -n p -i $indel -u $dup Evolved");

my $R = Statistics::R->new(); #bridge to use with R


###Run ClustalW



system ("clustalw -INFILE=Evolved -ALIGN -OUTFILE=1");
system ("clustalw -INFILE=1 -TREE");


#######Tree

my $cmds1 = <<EOF;
library("ape")
Standtree <- read.tree("StandardTREE")
QueTAL <- read.tree("Outputs/Output.tre")
#QueTAL2 <- read.tree("Outputs/Output2.tre")
Clust <- read.tree("1.ph")
Dist <- dist.topo(QueTAL,Standtree)
Dist2 <- dist.topo(Clust,Standtree)
EOF

my $run1 = $R->run($cmds1);
my $dist = $R->get('Dist');
my $dist2 = $R->get('Dist2');

push @arrayquetal, $dist;
push @arrayclust, $dist2;
#print $OUT1 "$s\t$dist\t$dist2\n";
}

}
my $median = median(@arrayquetal);
my $mean   = mean(@arrayquetal);
my $stdev  = stddev(@arrayquetal);
my $median2 = median(@arrayclust);
my $mean2   = mean(@arrayclust);
my $stdev2  = stddev(@arrayclust);

print $OUT3 "$probmut\t$probdel\t$mean\t$median\t$stdev\tDistal\n$probmut\t$probdel\t$mean2\t$median2\t$stdev2\tClustal\n";

}


}

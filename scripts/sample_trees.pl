#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);  #$Bin is the directory with the script
use lib "$Bin";        #add bin to the library path
use shared;
use shared_dups;
use Data::Dumper;

my $MUSCLE = "muscle";
my $MRBAYES = "mb";

my $MINLEN = 10;
my $BURNIN = 250;
my $NGEN = 10000;
my $SAMPLEFREQ = 10;

my $usage = "
$0 sequences atoms

Pre-computes phylogenetic tree samples for individual atoms.
";

my $sequences = shift or die $usage;
my $atoms = shift or die $usage;

# read FASTA file
open FA,"<$sequences" or die "Cannot open $sequences";
my %seqs;
while (my ($name,$seq) = read_fasta(\*FA)) {
    $name =~ s/^>//;
    $seqs{$name} = ${$seq};
}
close FA;

# read atoms and extract sequences
my %atoms; my %minlen;
open AT,"<$atoms" or die "Cannot open $atoms";
while (my $line = <AT>) {
    chomp $line;
    my @parts = split " ",$line;
    die "Wrong format" unless @parts == 6;
    push @{$atoms{$parts[2]}},\@parts;

    my $len = $parts[5]-$parts[4];
    if (exists $minlen{$parts[2]}) {
	$minlen{$parts[2]} = $len if $minlen{$parts[2]} > $len;
    } else {
	$minlen{$parts[2]} = $len;
    }

}
close AT;

my $tempdir = temp_dir_name($sequences);

foreach my $atomtype (keys %atoms) {
    print STDERR "Atom $atomtype\n";
    my $ntaxa = scalar(@{$atoms{$atomtype}});

    if ($minlen{$atomtype} <= $MINLEN || $ntaxa == 1) {
	print STDERR "... too short; writing star tree ($ntaxa atoms, length $minlen{$atomtype})\n";
	my @allatoms;
	foreach my $atom (@{$atoms{$atomtype}}) {
	    push @allatoms,$atom->[1].":0.00001"
	}
	print join("\t",$atomtype,0,"(".join(",",@allatoms).")root;"),"\n";
	open OUT,">$tempdir/$atomtype.aln.fa";
	foreach my $atom (@{$atoms{$atomtype}}) {
	    # extract sequence
	    my $atomseq = substr($seqs{$atom->[0]},$atom->[4],
				 $atom->[5]-$atom->[4]);
	    if ($atom->[3] < 0) {
		# reverse complement
		$atomseq = reverse($atomseq);
		$atomseq =~ tr/ACGTacgt/TGCAtgca/;
	    }

	    # put into fasta file
	    print OUT ">",$atom->[1],"\n",$atomseq,"\n";
	}
	close OUT;
    } else {
	## extract the sequences
	open OUT,">$tempdir/$atomtype.fa";
	foreach my $atom (@{$atoms{$atomtype}}) {
	    # extract sequence
	    my $atomseq = substr($seqs{$atom->[0]},$atom->[4],
				 $atom->[5]-$atom->[4]);
	    if ($atom->[3] < 0) {
		# reverse complement
		$atomseq = reverse($atomseq);
		$atomseq =~ tr/ACGTacgt/TGCAtgca/;
	    }

	    # put into fasta file
	    print OUT ">",$atom->[1],"\n",$atomseq,"\n";
	}
	close OUT;

	## sequence alignment
	my_run("cd $tempdir; $MUSCLE -in $atomtype.fa -out $atomtype.aln.fa >/dev/null")
	    unless (-r "$tempdir/$atomtype.aln.fa");

	## create nexus file for Mr. Bayes
	open IN,"<$tempdir/$atomtype.aln.fa";
	my %atomseqs; my $atomlen;
	while (my ($name,$seq) = read_fasta(\*IN)) {
	    $name =~ s/^>//;
	    $atomseqs{$name} = ${$seq};
	    $atomlen = length(${$seq});
	}
	close IN;

	## use every segment twice if the number is <4;
	open OUT,">$tempdir/$atomtype.nex";
	print OUT "#NEXUS\n";
	print OUT "begin data;\n";
	print OUT "dimensions ntax=",$ntaxa < 4 ? 2*$ntaxa: $ntaxa,
	" nchar=",$atomlen,";\n";
	print OUT "format datatype=dna gap=- interleave=yes;\n";
	print OUT "matrix\n";

	my $LINELEN = 50;
	for(my $i=0; $i < $atomlen; $i+=$LINELEN) {
	    foreach my $atomid (keys %atomseqs) {
		print OUT $atomid,"\t",
		substr($atomseqs{$atomid},$i,$LINELEN),"\n";
		print OUT "aux_$atomid\t",
		substr($atomseqs{$atomid},$i,$LINELEN),"\n" if ($ntaxa < 4)
	    }
	    print OUT "\n";
	}


	print OUT ";\nend;\n";
	print OUT "begin mrbayes;\n";
	print OUT "set autoclose=yes nowarn=yes;\n";
	print OUT "lset nst=2;\n";
	print OUT "mcmc nruns=1 ngen=$NGEN samplefreq=$SAMPLEFREQ;\n";
	print OUT "sumt burnin=$BURNIN;\n";
	print OUT "sump burnin=$BURNIN;\n";
	print OUT "end;\n";
	close OUT;

	## run Mr. Bayes
	my_run("cd $tempdir; $MRBAYES $atomtype.nex >/dev/null")
	    unless (-r "$tempdir/$atomtype.nex.p");

	## parse Mr. Bayes output
	open P,"<$tempdir/$atomtype.nex.p" or die;
	open T,"<$tempdir/$atomtype.nex.t" or die;
	my %translate;
	while (my $cmd = nexus_line(\*T)) {
	    if ($cmd =~ /translate (.*);/) {
		my @parts = split ",",$1;
		foreach my $part (@parts) {
		    my @map = split " ",$part;
		    $translate{$map[0]} = $map[1];
		}
	    } elsif ($cmd=~ /tree rep.(\d*) = (.*);/) {
		my ($sampleid,$strtree) = ($1,$2);
		next unless $sampleid >= $BURNIN*$SAMPLEFREQ;
		my $tree = read_tree_from_string($strtree.";");

		# rename nodes
		for(my $i=scalar(@$tree)-1; $i>=0; $i--) {
		    $tree->[$i]{'name'} = $translate{$tree->[$i]{'name'}} if
			defined $translate{$tree->[$i]{'name'}};
		    collapse_node($tree,$i) if
			$tree->[$i]{'name'} =~ /^aux_/ ||
			(node_is_leaf($tree->[$i]) &&
			 $tree->[$i]{'name'} =~ /^_aux_/);
		}

		my $newstrtree = write_tree_to_string($tree);

		# read corresponding sample from P file
		my @parts;
		while (my $line = <P>) {
		    next unless $line =~ /^\d/;
		    chomp $line;
		    @parts = split " ",$line;
		    last if ($parts[0] == $sampleid);
		}

		shift @parts;
		my $lnl = shift @parts;

		print join("\t",$atomtype,$lnl,$newstrtree,@parts),"\n";
	    }
	}

	close T;
	close P;
    }
}

sub nexus_line {
    my ($file) = @_;
    my $outline = "";
    while (my $line = <$file>) {
	chomp $line;
	$outline .= $line;
	return $outline if $line =~ /;/
    }
    return $outline;
}




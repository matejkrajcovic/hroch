package shared;

use strict;
use POSIX;
use DBI;
use Data::Dumper;

BEGIN {
    use Exporter   ();
    our (@ISA, @EXPORT, @EXPORT_OK);
    @ISA = qw(Exporter);
    # symbols to export by default
    @EXPORT = qw(read_fasta normalize_fasta_name write_fasta 
		 read_gtf parse_gtf_line write_gtf
		 parse_genepred_line write_genepred
		 my_run 
		 db_connect
	       );
}

############################
sub my_run
{
    my ($run, $die) = @_;
    if(!defined($die)) { $die = 1; }

    my $short = substr($run, 0, 20);

    print STDERR $run, "\n";
    my $res = system("bash", "-c", $run);
    if($res<0) {
        die "Error in program '$short...' '$!'";
    }
    if($? && $die) {
        my $exit  = $? >> 8;
        my $signal  = $? & 127;
        my $core = $? & 128;

        die "Error in program '$short...' "
            . "(exit: $exit, signal: $signal, dumped: $core)\n\n ";
    }
}


############################
sub read_fasta {
    # Read one fasta sequence from the fasta file
    # Return undef, if no sequence found; 
    #        name and reference to sequence otherwise.
    
    my ($input) = @_;
    
    # read name of fasta sequence
    my $name;
    while(1) {
	my $line = <$input>;
	
	# end of file
	return unless defined $line;

	# skip empty lines
	next if $line =~ /^\s*$/;

	# parse the name
	$line =~ s/\s+$//;
	if($line =~ /^>/) {
	    $name = $line; 
	    last;
	}
	else { die "Incorrect fasta file '$line'."; }
    }

    # read fasta sequence
    my $sequence = "";
    while(1) {
	my $file_pos = tell($input);
	my $line = <$input>;
	
	# end of file
	last unless defined $line; 

	# if there is a beginning of a new sequence
	if($line =~ /^>/) {
	    seek($input, $file_pos, 0);
	    last;
	}

	# remove all whitespaces from line and append it to the sequence
	$line =~ s/\s//g;
	$sequence .= $line;
    }
    
    return ($name, \$sequence);
}

############################
sub normalize_fasta_name
{
    my ($name) = @_;

    $name =~ s/\s+$//;          # remove trailing white space
    $name =~ s/^>\s*//;         # remove > at the beginning of line
    $name =~ s/^(\S+)\s.*$/$1/; # remove everything after first space

    #if it has gi, take only number
    if ( $name =~ /gi\|(.*)\|/ ) {
        $name = $1;
    }

    return $name;
}

############################
sub write_fasta {
    my ($file, $name, $seq) = @_;

    print $file $name, "\n";
    my $n = length($$seq);
    my $linelen = 60;
    my $i=0;
    while($i<$n) {
        print $file (substr($$seq, $i, $linelen), "\n");
        $i+=$linelen;
    }
}

############################
sub read_gtf
{
    my ($in, $name) = @_;

    my @rows;

    while(1) {
	my $file_pos = tell $in;   # remember where the line starts
	
	my $line = <$in>;
	last if !defined($line);

	next if $line =~ /^\#/;    #skip comments
	$line =~ s/\s+$//;         #strip trailing whitespace
	next if $line eq '';       #skip empty lines

	my ($seq, $rest) = split "\t", $line, 2;
	if($seq ne $name) {
	    seek $in, $file_pos, SEEK_SET;   #return to the start of line
	    last; 
	}
	else {                     # add one more line to sequence
	    push @rows, $line;
	}
    }
    return (\@rows);
}

########################
sub parse_gtf_line
{
    my ($line) = @_;

    return if $line =~ /^\#/;    #skip comments
    $line =~ s/\s+$//;
    return if $line eq '';       #skip empty lines

    my %res;

    @res{'seqname', 'source', 'feature', 'start', 'end', 'score',
        'strand', 'frame', 'attributes'} = split "\t", $line;

    die "Bad gtf line '$line'\n " if(!defined $res{'frame'});

    if(defined($res{attributes})) {
	my @attr = split ";", $res{"attributes"};
	foreach my $attrib (@attr) {
	    if($attrib =~ /^\s*(\S+)\s+"(.*)"\s*$/) {
		if(!exists($res{$1})) {
		    $res{$1} = $2;
		}
	    }
	}
    }
    return \%res;
}

############################
sub write_gtf
{
    my ($output, $rows) = @_;

    my @gtf_col_keys = ("seqname", "source", "feature", "start",
                        "end", "score", "strand", "frame");
    my %gtf_mandatory_cols;
    foreach my $col (@gtf_col_keys) { $gtf_mandatory_cols{$col} = 1; }

    foreach my $row (@$rows) {

        my $sep = "";
        foreach my $col (@gtf_col_keys) {
	    die "Bad row (missing $col)", Dumper($row), " "
		unless defined $row->{$col};
	    print $output $sep, $row->{$col};
	    $sep = "\t";
        }

	# [attributes]
        foreach my $col (sort keys %$row) {
            if(!($col eq "attributes")
               && !(exists $gtf_mandatory_cols{$col})) {
                printf $output "%s%s \"%s\";", $sep,
                $col, $row->{$col};
                $sep = " ";
            }
        }
	print $output "\n";
    }
}


our @genepred_cols=('name', 'chrom', 'strand', 'txStart', 'txEnd', 
		    'cdsStart', 'cdsEnd', 'exonCount', 
		    'exonStarts', 'exonEnds');

our @genepred_ext_cols = (@genepred_cols, 'id', 'name2', 
			  'cdsStartStat', 'cdsEndStat', 'exonFrames');

#############################
sub parse_genepred_line {

    my ($line) = @_;

    return if $line =~ /^\#/;    #skip comments
    $line =~ s/\s+$//;
    return if $line eq '';       #skip empty lines


    my @parts = split "\t", $line;

    my %rec;
    if(@parts == scalar(@genepred_cols)) {
	@rec{@genepred_cols} = @parts;
    }
    elsif(@parts == scalar(@genepred_ext_cols)) {
	@rec{@genepred_ext_cols} = @parts;
    }
    else {
	die "Wrong format '$line' ";
    }

    #split exon arrays and check their length
    my @exonStarts = split ',', $rec{'exonStarts'};
    my @exonEnds = split ',', $rec{'exonEnds'};

    my $n = $rec{'exonCount'};
    die unless $n>0;
    die unless scalar(@exonStarts)==$n 
      && scalar(@exonEnds)==$n;
    die unless $exonStarts[0]==$rec{'txStart'} && $exonEnds[$n-1]==$rec{'txEnd'};

    #store them as arrays in the record
    $rec{'exonStarts'} = \@exonStarts;
    $rec{'exonEnds'} = \@exonEnds;

    #handle optional exon frames
    my @exonFrames;
    if(exists $rec{'exonFrames'}) {
	@exonFrames = split ',',  $rec{'exonFrames'};
	die unless scalar(@exonFrames)==$n;
	$rec{'exonFrames'} = \@exonFrames;
    }
    return \%rec;
}

#############################
sub write_genepred {
    my ($file, $genepred) = @_;

    my %print = %$genepred;
    $print{'exonStarts'} = join(",", @{$print{'exonStarts'}});
    $print{'exonEnds'} = join(",", @{$print{'exonEnds'}});
    if(exists $print{'exonFrames'}) {
	$print{'exonFrames'} = join(",", @{$print{'exonFrames'}});
    }

    my $n = scalar keys %print;
    my @print_list;
    if($n == scalar(@genepred_cols)) {
	@print_list = @print{@genepred_cols};
    }
    elsif($n == scalar(@genepred_ext_cols)) {
	@print_list = @print{@genepred_ext_cols};
    }
    else {
	die "Wrong number of keys: $n";
    }

    print $file join("\t", @print_list), "\n";
}

#############################
sub db_connect {
    my ($database) = @_;

    #connect to hgsql db using ~/.hg.conf
    local *IN;
    open IN, "<", "/home/bbrejova/.hg.conf" or die;
    
    my ($host, $username,$password);
    while(my $line = <IN>) {
	$line =~ s/\s+$//;  #trim whitespace
	my ($left,$right) = split '=', $line;
	if($left eq 'db.host') { $host = $right; }
	elsif($left eq 'db.user') { $username = $right; }
	elsif($left eq 'db.password') { $password = $right; }
    }
    close IN;
    die unless defined $host && defined $username && defined $password;

    my $dbh = DBI->connect_cached("dbi:mysql:$database;$host",
			       $username,$password)
	or die "Cannot connect to the database $database@$host";
    return $dbh;
}




1;

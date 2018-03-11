package shared_dups;

use strict;
use POSIX;
use Data::Dumper;
use Bio::TreeIO;
use IO::String;

BEGIN {
    use Exporter   ();
    our (@ISA, @EXPORT, @EXPORT_OK);
    @ISA = qw(Exporter);
    # symbols to export by default
    @EXPORT = qw(read_history parse_history_line read_tree write_history read_tree_from_string 
                 write_tree write_tree_to_string node_is_leaf node_degree
                 collapse_node collapse_tree find_leftmost_leaf
                 find_node_by_name
                 read_atomic_sequence format_atoms parse_atoms temp_dir_name
                 sample_uniform_integer
                 read_config
	       );
}

############################
sub temp_dir_name {
    my ($query) = @_;
    # gets name of the fasta and creates a directory that
    # reflacts this name. This is then used to store
    # temporary files of exonhunter.

    my $res = "dupstemp_" . $query;
    $res =~ s/\W+/-/g;
    if(!-d $res) {
	mkdir($res) or die "Cannot create temporary directory $res";
	print STDERR "Created temporary directory $res\n";
    }

    return $res;
}


################################
sub format_atoms {
    my ($atoms) = @_;

    my @str_atoms;
    foreach my $atom (@$atoms) {
	my $str;
	if(!defined $atom) { $str = '-'; }
	else {
	    $str = $atom->{'type'} . "." . $atom->{'id'};
	    if($atom->{'strand'} == -1) {
		$str = '-' . $str;
	    }
	}
	push @str_atoms, $str;
    }
    return join(";", @str_atoms);
}

#####################
sub write_tree_to_string {
    my ($tree) = @_;
    my $string_fh = new IO::String;
    write_tree($string_fh,$tree);
    return ${$string_fh->string_ref};
}

#####################
sub write_tree {
    my ($file, $tree) = @_;

    rec_write_tree($file, $tree, 0);
    print $file ";";
}

#####################
sub rec_write_tree {
    my ($file, $tree, $node) = @_;

    my $time_str = "";
    if(defined $tree->[$node]{'length'} 
       && defined $tree->[$node]{'parent'}) { 
	$time_str = ":" . $tree->[$node]{'length'};
    }
    my $name_str = "";
    if(defined $tree->[$node]{'name'}) {
	$name_str = $tree->[$node]{'name'};
    }

    if(@{$tree->[$node]{'children'}}) {
	print $file "(";
	for(my $i=0; $i<@{$tree->[$node]{'children'}}; $i++) {
	    if($i>0) {
		print $file ",";
	    }
	    rec_write_tree($file, $tree, $tree->[$node]{'children'}[$i]);
	}
	print $file ")";
    }
    printf $file "%s%s", $name_str, $time_str;
}


######################
sub read_tree_from_string {
    my ($string) = @_;
    my $string_fh = new IO::String($string);
    return read_tree($string_fh);
}

#####################
sub read_tree {
    my ($file) = @_;

    #returns reference $tree with usage
    # $tree->[$node]{'name'/'parent'/'length'/'children'/'depth'} 
    #length is branch length to parent
    #nodes are ordered in DFS order
    #depth is sum of branch lengths to parent

    #use bioperl to parse newick format
    my $input = new Bio::TreeIO(-fh   => $file,
				-noclose => 1,
				-format => "newick");
    my $tree = $input->next_tree;
    die "Bad tree" unless defined $tree;
    
    # get our simple data structure with nodes 0,1,2...
    my @tree;
    my %name2num;
    my $num = 0;
    my $auxnum = 0;
    for my $node ( $tree->get_nodes()) {
	#we need all nodes to have unique names
	if ($node->id() eq "") {
	    $node->id(sprintf("_aux_%04d",$auxnum));
	    $auxnum++;
	}
	my $name = $node->id();
	die "non-unique name '$name'" if exists $name2num{$name};
	$name2num{$name} = $num;

	my $parent_node = $node->ancestor();
	my $parent;
	if(defined $parent_node) {
	    my $parent_name = $parent_node->id();
	    die unless defined $parent_name && exists $name2num{$parent_name};
	    $parent = $name2num{$parent_name};
	}
	my %rec = ('parent'=>$parent, 
		   'name'=>$name);
	if(defined $node->branch_length()) {
	    $rec{'length'}=$node->branch_length();
	}
	push @tree, \%rec;
	$num++;
	die unless $num==@tree;
    }

    #fill in child info, depth
    for(my $node=0; $node<@tree; $node++) {
	$tree[$node]{'children'} = [];  #empty array of children
	if($node>0) {
            #add node to parent's children
	    my $parent = $tree[$node]{'parent'};
	    die unless defined $parent;
	    push @{$tree[$parent]{'children'}}, $node; 
	    #add branch length to parent's depth
	    if(defined $tree[$node]{'length'} 
	       && defined $tree[$parent]{'depth'}) {
		$tree[$node]{'depth'} 
		= $tree[$parent]{'depth'}+$tree[$node]{'length'};
	    }
	}
	else {
	    $tree[$node]{'depth'} = 0;  #root is in depth 0
	}
    }

    return \@tree;

    #http://www.bioperl.org/wiki/HOWTO:Trees
    #http://doc.bioperl.org/releases/bioperl-current/bioperl-live/Bio/Tree/NodeI.html
    #http://doc.bioperl.org/releases/bioperl-current/bioperl-live/Bio/Tree/Tree.html
    #http://doc.bioperl.org/releases/bioperl-current/bioperl-live/Bio/TreeIO.html
}


######################
sub node_is_leaf {
    my ($node_rec) = @_;

    return @{$node_rec->{'children'}}==0;
}

######################
sub node_degree {
    my ($node_rec) = @_;

    my $degree = @{$node_rec->{'children'}};
    if(defined $node_rec->{'parent'}) { $degree++; }
    return $degree;
}

######################
sub collapse_node {
    my ($tree, $node) = @_;

    #Remove $node from $tree. It should not be the root node.
    #Its children become children of its parent, 
    # their branch lengths increase accordingly.
    #If the node is a leaf, it is simply deleted.
    
    my $parent = $tree->[$node]{'parent'};
    die unless defined $parent;
    
    #find $node in $parent's children
    my $index;
    for(my $i=0; $i<@{$tree->[$parent]{'children'}}; $i++) {
	if($tree->[$parent]{'children'}[$i]==$node) {
	    die if defined $index;
	    $index = $i;
	}
    }
    die unless defined $index;

    #replace $node in $parent's child list by $node's children
    splice @{$tree->[$parent]{'children'}}, $index, 1, 
           @{$tree->[$node]{'children'}};

    #link $parent as parent of $node's children
    #update their branch length
    my $add_length = $tree->[$node]{'length'};
    foreach my $child  (@{$tree->[$node]{'children'}}) {
	die unless $tree->[$child]{'parent'} == $node;
	$tree->[$child]{'parent'} = $parent;
	$tree->[$child]{'length'} += $add_length;
    }
						
    $tree->[$node] = {'deleted'=>1};  #empty record
}

######################
sub collapse_tree {
    my ($tree, $mode) = @_;

    # Collapse (delete) all internal vertices of degree 2
    # as well as root of degree 1
    # for $mode eq 'unrooted' also collapse root of degree 2
    # Keep root at node 0, leaves also keep their numbers.

    # Check and collapse root of degree 1
    if(node_degree($tree->[0])==1) {
	my $child = $tree->[0]{'children'}[0];
	if(!node_is_leaf($tree->[$child])) {
	    #root of degree 1, its child not a leaf - collapse the child
	    collapse_node($tree, $child);
	}
    }
    
    #check and collapse root of degree 2 if unrooted
    if($mode eq "unrooted") {
	while(node_degree($tree->[0])==2) {
	    my $child1 = $tree->[0]{'children'}[0];
	    my $child2 = $tree->[0]{'children'}[1];
	    if(!node_is_leaf($tree->[$child1])) {
		collapse_node($tree, $child1);
	    }
	    elsif(!node_is_leaf($tree->[$child2])) {
		collapse_node($tree, $child2);
	    }
	    else {
		last;  #both children are leaves, nothing to collapse
	    }
	}
    }

    #collapse other nodes of degree 2
    for(my $node = 1; $node < @$tree; $node++) {
	next if $tree->[$node]{'deleted'};
	die unless defined $tree->[$node]{'parent'};
	if(node_degree($tree->[$node])==2) {
	    collapse_node($tree, $node);
	}
    }
}

######################
sub find_leftmost_leaf {
    my ($tree) = @_;

    #fill in the node number of leftmost leaf for each node
    #starting from leaves
    
    for(my $node = @$tree-1; $node>=0; $node--) {
	next if $tree->[$node]{'deleted'};
	if(node_is_leaf($tree->[$node])) {
	    $tree->[$node]{'leftmost_leaf'} = $node;
	}
	else {
	    my $child = $tree->[$node]{'children'}[0];  #leftmost child
	    die unless defined $child 
		&& defined $tree->[$child]{'leftmost_leaf'};
	    $tree->[$node]{'leftmost_leaf'} = $tree->[$child]{'leftmost_leaf'};
	}
    } 
}

######################
sub find_node_by_name {
    my ($tree, $name) = @_;

    for(my $node = 0; $node<@$tree; $node++) {
	next if $tree->[$node]{'deleted'};
	if(exists $tree->[$node]{'name'} 
	   && $tree->[$node]{'name'} eq $name) {
	    return $node;
	}
    }
    return;
}

#####################
sub read_history {
    my ($file) = @_;
    
    # $history->[$num]{'source_species'/'target_species'/'time'/
    #                  'type'/'id'/'event_count'/'source_atoms'/'target_atoms'}
    # Atoms are records 'type', 'id', 'strand', 'species' 
    # where 'strand' is 1 or -1
    # Undefined atoms are deleted 

    my @history;
    while(my $line = <$file>) {
	my $rec = parse_history_line($line);
	if(defined $rec) {
	    push @history, $rec;
	}
    }

    return \@history;
}

############################
sub parse_history_line {
    my ($line) = @_;

    my @parts = split ' ', $line;	
    return if @parts==0 || $parts[0] eq "==="; 
    
    die unless @parts>=8;
    my $source_atoms = parse_atoms($parts[6], $parts[0]);
    my $target_atoms = parse_atoms($parts[7], $parts[1]);
    my %rec;
    @rec{'source_species','target_species',
	 'time','type','strand','event_count'} = @parts[0..5];
    @rec{'source_atoms','target_atoms'} = ($source_atoms, $target_atoms);
    $rec{'comments'} = join(" ", @parts[8..$#parts]);
    return \%rec;
}


############################
sub parse_atoms {
    my ($str, $species) = @_;
    my @parts = split ';', $str;
    my @atoms;
    foreach my $atom (@parts) {
	#deal with deleted atoms
	if($atom eq '-') {
	    push @atoms, undef;
	    next;
	}
	#normal atoms
	my %rec = ('strand' => 1);
	if(substr($atom,0,1) eq '-') {
	    $rec{'strand'} = -1;
	    $atom = substr($atom, 1); #remove leading -
	}
	my ($type, $id) = split '\.', $atom;
	die unless defined $id;
	$rec{'type'} = $type;
	$rec{'id'} = $id;
	if(defined $species) {
	    $rec{'species'} = $species;
	}
	push @atoms, \%rec;
    }
    return \@atoms;
}

#####################
sub write_history {
    my ($out, $history) = @_;

    foreach my $event (@$history) {
	printf $out "%s\t%s\t%.6g\t%s\t%d\t%d\t%s\t%s",
	$event->{'source_species'}, $event->{'target_species'},
	$event->{'time'}, $event->{'type'}, $event->{'strand'},
	$event->{'event_count'},
	format_atoms($event->{'source_atoms'}),
	format_atoms($event->{'target_atoms'});
	if(exists $event->{'comment'}) {
	    print $out "\t", $event->{'comment'};
	}
	print $out "\n";
    }
}

############################
sub read_atomic_sequence {
    my ($in) = @_;

    # $species2atoms->{$species}[$i]{'start'/'end'/'strand'/'species'/'id'/'length'}
    # straight from the file

    my $species2atoms;
    my %seen_id;
    while(my $line = <$in>) {
        my @parts = split ' ', $line;
        die "Bad line '$line'"unless @parts==6;
        my %rec;
        @rec{qw/species id type strand start end/} = @parts;
        $rec{'length'} = $rec{'end'}-$rec{'start'};
        die if exists $seen_id{$rec{'id'}};
        $seen_id{$rec{'id'}}=1;
        push @{$species2atoms->{$rec{'species'}}}, \%rec;
    }
    
    return $species2atoms;
}


######################
sub sample_uniform_integer {
    my ($max) = @_;
    #uniform integer in {0,1,...$max-1}

    my $r = POSIX::floor(rand($max));
    die unless 0<=$r && $r<$max;
    return $r;
}

######################
sub read_config {
    my ($options, $file, $required_keys) = @_;

    #Config file has lines containing keyword and value seperated by whitespace
    #Lines starting with # are comments.
    
    my %double;
    while(my $line = <$file>) {
	#skip commented lines
	next if $line=~ /^\#/;
	my @parts = split ' ', $line;
	#skip emtpy lines
	next if @parts==0;
	die unless @parts==2;
	$double{$parts[0]}++ if exists $options->{$parts[0]};
	$options->{$parts[0]} = $parts[1];
    }

    #warning for redefined keywords
    if(%double) {
	my $num = scalar(keys %double);
	warn "$num keywords occur more than once:\n'" 
	    . join("', '", sort keys %double) . "'\n ";
    }

    #die if required keys missing
    foreach my $key (@$required_keys) {
	die "Cannot find option $key" unless exists $options->{$key};
    }

    #warn if non-required keys present
    my %required_hash;
    @required_hash{@$required_keys} = ();
    my %extra;
    foreach my $key (keys %$options) {
	if(!exists $required_hash{$key}) {
	    $extra{$key}++;
	}
    }
    #warning for redefined keywords
    if(%extra) {
	my $num = scalar(keys %extra);
	warn "$num non-required keywords:\n'" 
	    . join("', '", sort keys %extra) . "'\n ";
    }
}


1;

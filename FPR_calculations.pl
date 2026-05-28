#! /usr/local/bin/perl

#####
#####
##
##  This script is to calculate precision, recall and F score of predicted GO annotations
##
##   Inputs:
##   -p  predicted annotations (tab delimited file)
##   -t  existing experimental annotations on the date the predictions were made
##   -r  reference experimental annotations generated after the date the predictions were made.
##   -n  GO do not annotate list. Terms determined not to be used for annotation by the GOC
##   -g  A lookup file with GO term to all its parent terms.
##   -o  output mapping file for each predicted annotation
##   STDOUT output FPR for each protein, as well as the average (at the end of the file)
###
#####
#####

# get command-line arguments
use Getopt::Std;
getopts('o:O:i:p:t:r:g:n:e:vVh') || &usage();
&usage() if ($opt_h);         # -h for help
$outFile = $opt_o if ($opt_o);    # -o for output files
$inFile = $opt_i if ($opt_i);     # -i for input file
$prediction = $opt_p if ($opt_p); # -p for the prediction file
$existing = $opt_t if ($opt_t);   # -t for existing annotation file.
$reference = $opt_r if ($opt_r);  # -r for new annotations file (as reference)
$do_not_annotate = $opt_n if ($opt_n); # -n for the list of do_not_annotate terms list
$go_parent = $opt_g if ($opt_g);  # -g for go parent file
$errFile = $opt_e if ($opt_e);    # -e for (e)rror file (redirect STDERR)
$verbose = 1 if ($opt_v);         # -v for (v)erbose (debug info to STDERR)
$verbose = 2 if ($opt_V);         # -V for (V)ery verbose (debug info STDERR)


###################################
# Working on GO parents
###################################

print STDERR "Working on GO child parent term lookup file.\n";
my %go_parents;   # a child to parent GO lookup
my %go_child;     # a parent to child GO lookup
my %go;           # a GO ID to term lookup

open (GP, $go_parent);
while (my $line=<GP>){
    chomp $line;
   
    my ($child, $parent, $relation, $foo)=split(/\t/, $line);
    
    my ($child_name, $child_id);
    
    if ($child=~/^(.+)\((GO\:\d+)\)$/){
        $child_name=$1;
        $child_id = $2;
    }
    
    my ($parent_name, $parent_id);
    if ($parent=~/^(.+)\((GO\:\d+)\)$/){
        $parent_name = $1;
        $parent_id = $2;
    }
##
    my $aspect;
    if ($foo=~/^m/){
        $aspect = 'mf';
    }elsif ($foo=~/^b/){
        $aspect = 'bp';
    }elsif ($foo=~/^c/){
        $aspect = 'cc';
    }
    $go_parents{$child_id}{$parent_id}=1;
    $go{$child_id}="$aspect\t$child_name";
    $go{$parent_id}="$aspect\t$parent_name";
    $go_child{$parent_id}{$child_id}=1;
}
close (GP);

############################################################################################
# Working on the do_not_annotate list
#  This is a list of GO terms determined by the GO Consortium not to be used for direct
#  annotation. They should be excluded from the comparison.
############################################################################################

print STDERR "Working on the do_not_annotate file.\n";
my %do_not_annotate;
open (DNA, $do_not_annotate);
while (my $line =<DNA>){
    chomp $line;
    my ($id, $name) = split(/\t/, $line);
    $do_not_annotate{$id}=1;
}
close (DNA);

#######################################################################################
#  Parse existing annotations
#   The file should contain the annotations to the parent terms.
#   If not, the script should be modified to use the goparent to add all the parents
#######################################################################################

print STDERR "Working on the existing annotation file.\n";
my %existing;
open (EX, $existing);
while (my $line=<EX>){
    chomp $line;
    my ($id, $go, $type) = split(/\t/, $line);
    next if ($go=~/GO:0005515|GO:0005488/);  # filter out binding and protein binding
    $existing{$id}{$go}=1;
}
close (EX);

################################################################################################
# Working on the reference (new annotations)
#  - read in the new annotation file (all annotations generated after the date of the prediction
#  - remove the redundant annotations that have parent-child relationships. Keep only the most
#    specific terms.
#  - remove existing terms. This could happen when the reference or EC was updated, so the new
#    creation date was assigned
#  - remove do_not_annotate terms
################################################################################################

print STDERR "Working on the new annotation file (reference annotations).\n";
my %new;
open (NEW, $reference);
while (my $line=<NEW>){
    chomp $line;
    my ($id, $go) = split(/\t/, $line);
    if (exists $existing{$id}{$go}){   # exclude annotations already exist to the gene.
        print STDERR "\-\-$id $go already in existing_annotation, removed from reference.\n";
        next;
    }
        
    next if ($go=~/GO:0005515|GO:0005488/);  # filter out binding and protein binding
    if (exists $do_not_annotate{$id}{$go}){  #remove do_not_annotate terms
        print STDERR "\-\-$id $go in the do_not_annotate list, removed from reference.\n";
        next;
    }
    $new{$id}{$go}=1;
}
close (NEW);

# Remove redundant terms.
# If a gene is annotated to terms that have parent child relationship, only the most specific
#  term should to retained.

my %reference = &removeRedundancy(\%new, \%go_parents, 'reference');

################################################################################################
# Working on the prediction file
#  - read in the "existing annotation" file that contains exp annotations on the date the
#    predictions were made
#  - remove annotations already in the "existing annotations" file
#  - remove do_not_annotate terms
#  - remove the redundant annotations that have parent-child relationships. Keep only the most
#    specific terms.
################################################################################################

print STDERR "Working on the predicted annotation file.\n";
my %hash;
open (FH, $prediction);
while (my $line=<FH>){
    chomp $line;
    my ($id, $go, $score) = split(/\t/, $line);
    next if ($go=~/GO:0005515|GO:0005488/);  # filter out binding and protein binding
    next if ($parent_go=~/GO:0008150|GO:0005575|GO:0003674/); # filter out the root terms
    if (exists $existing{$id}{$go}){   # exclude annotations already exist to the gene.
        print STDERR "\-\-$id $go already in existing_annotation, removed from predicted.\n";
        next;
    }
    if (exists $do_not_annotate{$id}{$go}){    #remove do_not_annotate terms
        print STDERR "\-\-$id $go in the do_not_annotate list, removed from predicted.\n";
        next;
    }
    $hash{$id}{$go}=1;
}
close (FH);

# Remove redundant annotations

my %predicted = &removeRedundancy(\%hash, \%go_parents, 'predicted');


##################################################################
# Mapping annotations between those the predicted and the reference
##################################################################

print STDERR "Mapping annotations.\n";
my %prediction_map;
my %prediction_map_count;  # This is to count the mapped terms.
my %prediction_nomap;  # No new exp for the gene. Not mapped.
my %reference_map;
my %reference_nomap;

foreach my $id (keys %predicted){
    if (exists $reference{$id}){
        foreach my $go_p (keys %{$predicted{$id}}){
            my %hash;
            foreach my $go_r (keys %{$reference{$id}}){
                #  my $type = $reference{$id}{$go_r};
                if ($go_p eq $go_r){
                    $hash{'direct'}{$go_r}=1;   # when the go terms in the predicted and reference are identical.
                }elsif (exists $go_parents{$go_p}{$go_r}){
                    $hash{'related'}{$go_r}=1;  # when the go term in the predicted is more specific to the one in the reference.
                }elsif (exists $go_parents{$go_r}{$go_p}){
                    $hash{'true'}{$go_r}=1;     # when the go term in the predicted is more general to the one in the reference.
                }
            }
            
            my $map;
            if (exists $hash{'direct'}){
                my @go_rs = keys %{$hash{'direct'}};
                my $go_rs = join("\;", @go_rs);
                $prediction_map{$id}{'direct'}{$go_p}=$go_rs;
                foreach my $go_r (@go_rs){
                    $reference_map{$id}{$go_r}=1;
                }
            }elsif (exists $hash{'true'}){
                my @go_rs = keys %{$hash{'true'}};
                my $go_rs = join("\;", @go_rs);
                $prediction_map{$id}{'true'}{$go_p}=$go_rs;
                $prediction_map_count{$id}{'true'}{$go_rs}=1; # when multiple more general GO terms are mapped to a specific exp GO term, they will only treated as on prediction.
                foreach my $go_r (@go_rs){
                    $reference_map{$id}{$go_r}=1;
                }
            }elsif (exists $hash{'related'}){
                my @go_rs = keys %{$hash{'related'}};
                my $go_rs = join("\;", @go_rs);
                $prediction_map{$id}{'related'}{$go_p}=$go_rs;
                foreach my $go_r (@go_rs){
                    $reference_map{$id}{$go_r}=1;
                }
            }else{
                $prediction_map{$id}{'unrelated'}{$go_p}=1;  # The go terms in the reference are not related to the ones in the predicted.
            }
        }
                
    }else{
        foreach my $go_p (keys %{$predicted{$id}}){
            $prediction_nomap{$id}{$go_p}=1;
        }
    }
}

foreach my $id (keys %reference){
    foreach my $go (keys %{$reference{$id}}){
        next if (exists $reference_map{$id}{$go});
        $reference_nomap{$id}{$go}=1;
    }
}
 
###############################################################
#@ Calculate the precision, recall and F score
###############################################################

print STDERR "Calculating precision, recall and F scores.\n";
my %precision;
my %recall;
foreach my $id (keys %prediction_map){
    my $E = keys (%{$prediction_map{$id}{'direct'}}); # Number of exactly matching annotations between the reference and the predicted set.
    my $L = keys (%{$prediction_map_count{$id}{'true'}}); # Number of predicted annotations that are less specific than an experimental annotation
    my $M = keys (%{$prediction_map{$id}{'related'}}); # Number of predicted annotations that are more specific than an experimental annotation
    my $A = keys (%{$prediction_map{$id}{'unrelated'}}); # Number of predicted annotations that are not related to any experimental annotations for the protein
    my $Z; #Number of experimental annotations that are not related to any predicted annotations for the protein
    if (exists $reference_nomap{$id}){
        $Z = keys (%{$reference_nomap{$id}});
    }
    
    my $precision =  ($E + 0.75*$L + 0.5*$M)/($E + $L + $M + $A);
    my $recall = ($E + 0.75*$L + 0.5*$M)/($E + $L + $M + $Z);
    
    $precision{$id}=$precision;
    $recall{$id}=$recall;
    
    $sum_precision+=$precision;
    $sum_recall+=$recall;
}

my $n_gene_mapped_predictions = keys (%{prediction_map});  # total number of proteins with at least one predicted GO term and at least one experimental GO term
my $n_gene_reference = keys (%reference);  # total number of proteins with at least one experimental GO term

my $ave_precision = ($sum_precision/$n_gene_mapped_predictions);
my $ave_recall = ($sum_recall/$n_gene_reference);

my $f_score = (2*($ave_precision)*($ave_recall)/($ave_precision + $ave_recall));

####################################
# print out results
####################################

open (OUT, ">$outFile");
print OUT "id\tpredicted\tmap type\treference\n";
foreach my $id (keys %prediction_map){
    foreach my $type (keys %{$prediction_map{$id}}){
        foreach my $go_p (keys %{$prediction_map{$id}{$type}}){
            my $go_r = $prediction_map{$id}{$type}{$go_p};
            $go_r = "" unless ($go_r=~/^GO/);
            print OUT "$id\t$go_p\t$type\t$go_r\n";
        }
    }
}

foreach my $id (keys %prediction_nomap){
    foreach my $go (keys %{$prediction_nomap{$id}}){
#        print OUT "$id\t$go\tno_map\t\n";
    }
}
close (OUT);

print "Average precision\t$ave_precision\n";
print "Average recall\t$ave_recall\n";
print "F_score\t$f_score\n";

print "\n\n";
print "Precision and Recall for individual proteins\n\n";
print "id\tpredicion\trecall\n";
foreach my $id (keys %precision){
    my $precision = $precision{$id};
    my $recall = $recall{$id};
    print "$id\t$precision\t$recall\n";
}



##############################
# Subroutines
##############################

sub removeRedundancy{
    my ($href, $go_parents_href, $type) = @_;
    
    my %hash;
    foreach my $id (keys %{$href}){
        my @gos = keys %{$href->{$id}};
        
        my @specific_gos;
        foreach my $go (@gos){
            my $is_redundant = 0;
            for my $other (@gos) {
                next if $other eq $go;
                # If $go appears in the parents of $other, it is a parent/ancestor -> redundant
                if (exists $go_parents_href->{$other}->{$go}) {
                    $is_redundant = 1;
                    print STDERR "\-\-$id $go is a parent GO term to $other, removed from $type.\n";
                    last;
                }
            }
            push @specific_gos, $go unless $is_redundant;
            
        }
        foreach my $go_s (@specific_gos){
            $hash{$id}{$go_s}=1;
        }
    }
    return %hash;
}





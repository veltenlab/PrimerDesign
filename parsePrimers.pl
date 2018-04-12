use strict;
use warnings;
#Aim: Parse primer3 output into an R-friendly format

my $id;
my @left;
my @right;

while (my $line = <>) {
    if ($line =~ /SEQUENCE_ID=(.+)/) {
       if ($id) {
        print $id . ";" . join(",",@left) . ";" . join(",",@right) . "\n"; 
       } 
       @left = (); @right = (); $id = $1;
    }
    if ($line =~ /PRIMER_LEFT_.+_SEQUENCE=(.+)/) {
        push(@left,$1);
    }
    if ($line =~ /PRIMER_RIGHT_.+_SEQUENCE=(.+)/) {
        push(@right,$1);
    }
}
print $id . ";" . join(",",@left) . ";" . join(",",@right) . "\n";

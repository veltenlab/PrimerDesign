use strict;
use warnings;
#Aim: Parse primer3 output into an R-friendly format

my ($id,$complany,$complend);

while (my $line = <>) {
    if ($line =~ /SEQUENCE_ID=(.+)/) {
       if ($id) {
        print "$id,$complany,$complend\n";
       } 
       $complany = ""; $complend=""; $id = $1;
    }
    if ($line =~ /PRIMER_PAIR_.+_COMPL_ANY_TH=(.+)/) {
        $complany = $1;
    }
    if ($line =~ /PRIMER_PAIR_.+_COMPL_END_TH=(.+)/) {
         $complend = $1;
    }
}
print "$id,$complany,$complend\n";

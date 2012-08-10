#!/usr/bin/env perl
use warnings;
use strict;
#
# Albert Danial     January 23 2001  
#
my ($usage, $script, );
($script = $0) =~ s{.*/}{};  # extract script's basename
$usage = "
Usage: $script  <t|b>  <dn|s1|s2>
               Arg 1:
                   t  = text file
                   b  = binary file

               Arg 2:
                   dn = dense matrix
                   s1 = sparse, BIGMAT = FALSE
                   s2 = sparse, BIGMAT = TRUE

";
die $usage unless @ARGV == 2;
my $file_format = shift @ARGV;
my $mat_format  = shift @ARGV;
die "first  arg must be t or b\n"       unless $file_format =~ /^[tb]$/;
die "second arg must be dn, s1 or s2\n" unless $mat_format  =~ /^(dn|s1|s2)$/;

my @template = <DATA>;
my $T = join("", @template);

$T =~ "ASSIGN output4='mat_${file_format}_${mat_format}.op4',";
$T =~ s/>FFORM</$file_format/mg;
$T =~ s/>MFORM</$mat_format/mg;
if ($file_format eq "t") {        # text   file
    $T =~ s/>UFORM<//mg;
} else {                          # binary file
    $T =~ s/>UFORM</UN/mg;
}
if      ($mat_format eq "dn") {   # dense
    $T =~ s/>UNIT_SIGN<//mg;
    $T =~ s/>BIGMAT</FALSE/mg;
} elsif ($mat_format eq "s1") {   # sparse, BIGMAT=FALSE
    $T =~ s/>UNIT_SIGN</-/mg;
    $T =~ s/>BIGMAT</FALSE/mg;
} elsif ($mat_format eq "s2") {   # sparse, BIGMAT=TRUE
    $T =~ s/>UNIT_SIGN</-/mg;
    $T =~ s/>BIGMAT</TRUE/mg;
}

my $outfile = "mk_op4.dat";
open(OUT, ">$outfile") or die "Cannot write to $outfile $!\n";
print OUT $T;
close(OUT);

print "Wrote $outfile.\n";

__END__
$ two steps:  1. uncomment one of the four ASSIGN statements
$             2. set IUnit =-12 for sparse,  IUnit = 12 for dense
$
$ ASSIGN output4='mat_ascii_sparse.op4', UNIT=12,FORMATTED,DELETE $ ASCII
$ ASSIGN output4='mat_ascii_dense.op4', UNIT=12,FORMATTED,DELETE $ ASCII
$ ASSIGN output4='mat_binary_sparse.op4', UNIT=12,UNFORMATTED,DELETE $ binary
$ ASSIGN output4='mat_bin_dense_i386.op4', UNIT=12,UNFORMATTED,DELETE $ binary
ASSIGN output4='mat_>FFORM<_>MFORM<.op4', UNIT=12,>UFORM<FORMATTED,DELETE
DIAG 8,56
TIME 3
SOL 100
compile userdmap souin=mscsou list noref $
alter 1,2 $
TYPE PARM,, I,N, ITape, IUnit, Digits, iPrec, Seed $ integer
TYPE PARM,, I,N, nRows, nCols, P5, P6, P7, P8, P9  $ integer
TYPE PARM,,LOGICAL,N,BigMat                        $ logical
TYPE PARM,,CS,N, Imagn, MoneI                      $ complex, single
$
IUnit  =>UNIT_SIGN<12      $ unit number, >0 = dense   <0 = sparse
BigMat =>BIGMAT<   $ only applies if sparse
$
$ Make some matrices
$
matgen /Eye10/1/10                    $ [I],         10 x 10
matgen /Eye5/1/  5                    $ [I],          5 x  5
matgen /Low/4/7/5/ / / /2/1/5         $ [pattern 1],  5 x  7
Seed   = 314159
Imagn  = CMPLX( 0.0, 1.0)             $  0 + i
MoneI  = CMPLX(-1.0, 1.0)             $ -1 + i
iPrec  = 1
matgen /Rnd1RS/5/4/4/iPrec/Seed////   $ [random  1],  4 x  4,  RS
iPrec  = 2
matgen /Rnd1RD/5/4/4/iPrec/Seed////   $ [random  2],  4 x  4,  RD
iPrec  = 1
add    Rnd1RS,/Rnd1CS/MoneI///        $ [random  3],  4 x  4,  CS
iPrec  = 2
add    Rnd1RD,/Rnd1CD/MoneI///        $ [random  4],  4 x  4,  CD
iPrec  = 2
matgen /Null/7/3/3/0/iPrec////        $ [0],  3 x  3,  RD
add    Eye5,/Eye5CD/MoneI///          $ [I - iI],  5 x  5,  CD
iPrec  = 2
iPrec  = 1
nRows  = 20
nCols  = 30
P5     =  3
P6     =  4
P7     =  3
P8     = 12
P9     =  2
matgen /Strings/4/nRows/nCols/iPrec/P5/P6/P7/P8/P9/  $ [pattern]
$
$ Write them to .op4
$
ITape  = 0       $ 0=append, -1=rewind before write, -2=end file and
                 $      rewind after write,  -3=both
Digits = 16       $ examples: 9 -> 1PE16.9 (default),   16 -> 1PE23.16
output4 Eye10,Low,,,//ITape/IUnit//BigMat/Digits $  
output4 Rnd1RS,Rnd1RD,Rnd1CS,Rnd1CD,// ITape/IUnit//BigMat/Digits $
Digits = 15       $ examples: 9 -> 1PE16.9 (default),   16 -> 1PE23.16
output4 Null,Strings,Eye5CD,,// ITape/IUnit//BigMat/Digits $
$ matprn Eye10,Low// $
endalter
cend
begin bulk
enddata

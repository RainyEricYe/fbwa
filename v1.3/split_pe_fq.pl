#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;
use FindBin qw/$Bin/;
use lib "/share/backup/yerui/PerlModule/my_pm";
use Iopen;

my ( $Help, $Line, );
GetOptions(
    'help|?' => \$Help,
    'line=i' => \$Line,
);
die `pod2text $0` if $Help or @ARGV < 2 or !$Line;

my ($f1, $f2, $p1, $p2) = @ARGV;

my $outdir = dirname($p1);
system "mkdir -p $outdir" unless -d $outdir;
$outdir = abs_path($outdir);

$Line *= 4; # 4 lines for a reads

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my ($IN1, $IN2);
iopen(\$IN1, "<", $f1);
iopen(\$IN2, "<", $f2);

my $part = 0;
my ($OUT1, $OUT2, $of1, $of2);
add_p();

my $n = 0;
my $cnt = 0;
while (1) {
    my $r1 = <$IN1>;
    my $r2 = <$IN2>;

    last if !$r1 or !$r2;

    ++$n;
    ++$cnt;
    print $OUT1 "$r1";
    print $OUT2 "$r2";

    if ( $cnt == $Line ) {
        $cnt = 0;
        add_p();
    }
}

map {close $$_} (\$IN1, \$IN2, \$OUT1, \$OUT2);

if ( $cnt == 0 ) {
    system "rm -f $of1 $of2";
    --$part;
}

print "total $n lines, $part parts\n";

&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub add_p {

    ++$part;
    $of1 = "$p1.$part.fq.gz";
    $of2 = "$p2.$part.fq.gz";

    iopen(\$OUT1, ">", $of1);
    iopen(\$OUT2, ">", $of2);
}

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900, $t[4] + 1, @t[3,2,1,0], $_[0];
}

__END__

=head1 Function

=head1 Usage
    perl $0 in1.fq(.gz) in2.fq(.gz) out1.prefix out2.prefix

=head1 Options
    -h             help
    -line  [int]   [int] reads as a part to ouput

=head1 Author
    yerui;    yerui@genomics.cn

=head1 Version
    v1.0;    2012-5-21

=cut

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;
use FindBin qw/$Bin/;

my ( $Help, $Pair, $Clean, $Line, $Verbose, $Queue, $Check, $Maxjob, ) =
   ( 0,     1,     0,      2**19, 0, " -q st.q -P st_cancer ",  0,  500, );

my $FQ_MEM    = "800M";
my $ALN_MEM   = "4G";
my $SAMPE_MEM = "4G";
my $MERGE_MEM = "4G";

my $Aln_opt = "-R 15 -I -L -m 1000000 -e 50";
my $Ref = "/ifs1/ST_CANCER/SHARE/Database/HumanHg19/Reference/hg19_new/hg19.fa";

GetOptions(
    'help'      => \$Help,
    'pair'      => \$Pair,
    'clean'     => \$Clean,
    'line=i'    => \$Line,
    'aln_opt=s' => \$Aln_opt,
    'ref=s'     => \$Ref,
    'queue=s'   => \$Queue,
    'maxjob=i'  => \$Maxjob,
    'verbose'   => \$Verbose,

    'fq_mem=s'    => \$FQ_MEM,
    'aln_mem=s'   => \$ALN_MEM,
    'sampe_mem=s' => \$SAMPE_MEM,
    'merge_mem=s' => \$MERGE_MEM,

    'check' => \$Check,
);

die `pod2text $0` if $Help or @ARGV < 2;

vshow("use reference $Ref");

my ( $fq_f_list, $outdir ) = @ARGV;

mkdir $outdir unless -d $outdir;
$outdir = abs_path($outdir);

my $split_dir = "$outdir/split_fq";
my $bam_dir   = "$outdir/bam";
my $sh_dir    = "$outdir/sh";

make_dir( $split_dir, $bam_dir, $sh_dir );

my $split_prog = "$Bin/split_pe_fq.pl";
my $filt_prog  = "$Bin/filt_adapter";
my $bwa        = "$Bin/bwa";
my $samtools   = "$Bin/samtools";
my $qq         = "$Bin/qinfo.pl";

check_prog( $split_prog, $filt_prog, $bwa, $samtools, $qq );

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

&showLog("START");

my $qsubed_job = 0;    # used in mqsub. qsub one job, $qsubed_job++
my %jobs;              # all jobs id in queue.

if ($Pair) {
    my @fq_fs = ();

    read_fq_flist( $fq_f_list, \@fq_fs );
    fq_is_pair( \@fq_fs );

    if ($Check) {
        kill_aln_jobs($sh_dir);
    }
    else {
        if ($Clean) {    # input clean fq files
            vshow("split pair end fq");
            split_pe_fq_files( \@fq_fs );
        }
        else {
            vshow("filt and split fq");
            filt_adapter_and_split_pe_fq_files( \@fq_fs );
        }
    }

    vshow("start align");

    bwa();

    if ($Check) {
        if ( $qsubed_job == 0 ) {
            showLog("All jobs finished correctly.");
        }
        else {
            showLog("Total $qsubed_job jobs have been qsubed this time.");
        }
    }
    else {
        showLog("Total $qsubed_job have been qsubed.");
    }
}
else {
    warn "under developing ...";
}

&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# kill all jobs in $outdir/sh, except *split_fq.sh
sub kill_aln_jobs {
    my ( $d, ) = @_;

    my $user = `whoami`;
    chomp $user;

    for my $j (`qstat -u $user | $qq -id | grep $d | grep -v split_fq`) {
        chomp $j;

        my ( $id, $sh ) = split /\s+/, $j;

        system "rm -f $sh.* && qdel $id";
        vshow("rm -f $sh.* && qdel $id");
    }

}

sub make_dir {
    for my $d (@_) {
        mkdir $d unless -d $d;
    }
}

sub check_prog {
    for my $p (@_) {
        die "cannot execute $p" if !-x $p;
    }
}

sub read_fq_flist {
    my ( $f, $r ) = @_;

    open my $LIST, $f or die $!;
    while ( my $k = <$LIST> ) {
        chomp $k;

        if ( !-f $k or !-s $k or $k !~ /\.fq(\.gz)?$/ ) {
            die "fq file error: $k";
        }

        push @$r, $k;
    }
}

sub fq_is_pair {
    my ( $r, ) = @_;
    my @t = @$r;

    while (@t) {
        my $f1 = shift @t;
        my $f2 = shift @t;

        if ( !$f2 ) {
            die "not pair $f1";
        }

        my $lane1 = ( split /_1/, basename($f1) )[0];
        my $lane2 = ( split /_2/, basename($f2) )[0];

        if ( !$lane1 or $lane1 ne $lane2 ) {
            die "not pair $f1";
        }
    }

    return 1;
}

sub split_pe_fq_files {
    my ( $r, ) = @_;

    while (@$r) {
        my $f1 = shift @$r;
        my $f2 = shift @$r;

        my $lane    = ( split /\//, $f1 )[-2];
        my $outd    = "$split_dir/$lane";
        my $fq_sh_d = "$sh_dir/$lane";

        make_dir( $outd, $fq_sh_d );

        open my $SH, ">$fq_sh_d/split_fq.sh" or die $!;
        print $SH
"time $split_prog -l $Line $f1 $f2 $outd/a $outd/b > $outd/stat && echo job-done\n";
        close $SH;

        mqsub( "$fq_sh_d/split_fq.sh", $FQ_MEM );
    }
}

sub filt_adapter_and_split_pe_fq_files {

    my ( $r, ) = @_;

    while (@$r) {
        my $f1 = shift @$r;
        my $f2 = shift @$r;

        my $lane    = ( split /\//, $f1 )[-2];
        my $outd    = "$split_dir/$lane";
        my $fq_sh_d = "$sh_dir/$lane";

        make_dir( $outd, $fq_sh_d );

        my $dir = dirname($f1);
        my @ads = glob "$dir/*adapter.list*";
        chomp @ads;

        if ( @ads != 2 ) {
            die "ERROR num of adapter list in $dir";
        }

        open my $SH, ">$fq_sh_d/filt.split_fq.sh" or die $!;
        print $SH
"time $filt_prog -p $Line $f1 $f2 @ads $outd/a $outd/b > $outd/stat && echo job-done\n";
        close $SH;

        mqsub( "$fq_sh_d/filt.split_fq.sh", $FQ_MEM );
    }

}

sub bwa {

    my %hid;    # hold jid for each lane_dir

    my @dirs = glob "$split_dir/*/";
    chop @dirs;

    my %h;      # key => dirs, value => which part has been done
    map { $h{$_} = 1 } @dirs;

    while (1) {
        my $flag = 1;    # whether all lane finished

      LANE:
        for ( my $i = 0 ; $i < @dirs ; $i++ ) {
            my $d = $dirs[$i];

            if ($d) {
                $flag = 0;
            }
            else {
                next LANE;
            }

            while ( new_fq_ok( $d, $h{$d} ) ) {
                my $sam_id = aln_sampe_sort( $d, $h{$d} );

                $hid{$d} .= ( $sam_id ? "$sam_id," : "" );
                $h{$d}++;
            }

            if ( split_done($d) ) {
                merge_rmdup( $d, $h{$d} - 1, $hid{$d} );
                $dirs[$i] = 0;    # finish split & aln
            }
        }

        last if $flag;

        sleep 10;
    }
}

# param split_dir/lane
# param part_num of fq
sub new_fq_ok {
    my ( $d, $n ) = @_;

    my $a = "$d/a." . ( $n + 1 ) . ".fq.gz";
    my $b = "$d/b." . ( $n + 1 ) . ".fq.gz";

    if (   ( -f $a && -f $b )
        || ( split_done($d) && -f "$d/a.$n.fq.gz" && -f "$d/b.$n.fq.gz" ) )
    {
        return 1;
    }
    else {
        return 0;
    }
}

# param split_dir/lane
# param part_num of fq
# return jid of sampe

sub aln_sampe_sort {

    my ( $d, $n ) = @_;

    my $sh_d = "$sh_dir/" . basename($d);
    make_dir($sh_d);

    my $aln_d   = "$sh_d/aln";
    my $sampe_d = "$sh_d/sampe";
    make_dir( $aln_d, $sampe_d );

    my $sampe_sh = "$sampe_d/s$n.sh";

    vshow("aln sampe $d");

    for my $t (qw/a b/) {
        open my $SH, ">$aln_d/$t.$n.sh" or die $!;
        print $SH "time $bwa aln " . $Aln_opt
          . " $Ref $d/$t.$n.fq.gz -f $d/$t.$n.fq.gz.sai && echo job-done\n";
        close $SH;
    }

    open my $SH, ">$sampe_sh" or die $!;
    print $SH
"time ($bwa sampe $Ref $d/a.$n.fq.gz.sai $d/b.$n.fq.gz.sai $d/a.$n.fq.gz $d/b.$n.fq.gz | $samtools view -S - -bo $d/$n.bam && $samtools sort -m 3500000000 $d/$n.bam $d/$n.srt && rm -f $d/$n.bam && echo job-done)\n";
    close $SH;

    my @ids = ();

    my $a_id = reqsub( "$aln_d/a.$n.sh", $ALN_MEM );
    push @ids, $a_id if $a_id;

    my $b_id = reqsub( "$aln_d/b.$n.sh", $ALN_MEM );
    push @ids, $b_id if $b_id;

    return reqsub( $sampe_sh, $SAMPE_MEM, ( join ",", @ids ) );
}

# param split_dir/lane
sub split_done {
    my ( $d, ) = @_;

    return 0 if !-f "$d/stat";

    my $n = ( split /\s+/, `wc -l $d/stat` )[0];

    if ($Clean) {
        $n ? return 1 : 0;
    }
    else {
        $n > 1 ? return 1 : 0;
    }
}

# param split_dir/lane
# param num of parts
sub merge_rmdup {
    my ( $d, $n, $hid ) = @_;

    my @bam_f;
    map { push @bam_f, "$d/$_.srt.bam" } ( 1 .. $n );

    my $lane = basename($d);

    my $outd = "$bam_dir/$lane";
    my $sh   = "$sh_dir/$lane/merge_rmdup.sh";

    make_dir( $outd, "$sh_dir/$lane" );

    open my $SH, ">$sh" or die $!;

    if ( @bam_f > 1 ) {
        print $SH
"time ($samtools merge -f $outd/srt.bam @bam_f && $samtools rmdup $outd/srt.bam $outd/rmd.bam && rm -f $outd/srt.bam && $samtools index $outd/rmd.bam && echo job-done)";
    }
    else {
        print $SH
#"time (mv -f $bam_f[0] $outd && $samtools rmdup $outd/srt.bam $outd/rmd.bam && rm -f $outd/srt.bam && $samtools index $outd/rmd.bam && echo job-done)";
        "time ($samtools rmdup $bam_f[0] $outd/rmd.bam && $samtools index $outd/rmd.bam && echo job-done)"; # thank zhouquau
    }

    close $SH;

    chop $hid if $hid =~ /,$/;

    return reqsub( $sh, $MERGE_MEM, $hid );
}

sub mqsub {
    my ( $sh, $mem, $hid ) = @_;

    my $id;
    my ( $nsh, $dir ) = fileparse($sh);
    chdir $dir;

    while (1) {
        my $cmd =
            "qsub -cwd -l vf=$mem $Queue "
          . ( $hid   ? " -hold_jid $hid " : " " ) . "$nsh";
        my $info = `$cmd`;
        $id = ( $info =~ /Your job (\d+) / )[0];
        if ( defined $id ) {
            print "$id --> $cmd\n";
            last;
        }

        sleep 60;
    }

    $qsubed_job++;

    $jobs{$id} = 1;
    while ( scalar( keys %jobs ) >= $Maxjob ) {
        sleep 120;

        for my $j ( keys %jobs ) {
            if ( !`qstat -j $j` ) {
                delete $jobs{$j};
            }
        }
    }

    return $id;
}

sub reqsub {
    my ( $sh, $mem, $hid ) = @_;

    my @log = glob "$sh.*";
    if ( @log == 2 && `grep job-done $log[1]` ) {
        return 0;
    }
    else {
        if (@log) {
            system "rm -f @log";
            vshow("reqsub rm -f @log");
        }

        return mqsub( $sh, $mem, $hid );
    }
}

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900,
      $t[4] + 1, @t[ 3, 2, 1, 0 ], $_[0];
}

sub vshow {
    showLog( $_[0] ) if $Verbose;
}

__END__

=head1 Function
    Input raw fq file list, output bam results by lane.
    This program split fq files into small pieces to do alignment.

=head1 Usage
    perl $0 $fq_f_list $outdir

=head1 Options
    -help          help
    -pair          pair end sequence [true]
    -clean         fq files are cleaned (do not need filt adapter) [false]
    -line    [i]   put [int] reads as a part [2**19]
    -aln_opt [s]   option for bwa_aln. [ -R 15 -I -L -m 1000000 -e 50 ]
    -ref     [s]   reference.fa   ["/ifs1/ST_CANCER/SHARE/Database/HumanHg19/Reference/hg19_new/hg19.fa"]
    -queue   [s]   queue request [" -q st.q -P st_cancer "]
    -verbose       show verbose log
    -maxjob  [i]   maximum num of jobs in queue. [500]
    -check         if some job failed, except *split_fq.sh, you can restart the main perl with this option.
                   You can check the result with this option.

=head1
    -fq_mem    [s]      memory request for filt.split_fq.sh or split_fq.sh [800M]
    -aln_mem   [s]      memory request for bwa aln [4G]
    -sampe_mem [s]      memory request for bwa sampe [4G]
    -merge_mem [s]      memory request for samtools merge [4G]

=head1 Author
    yerui;    yerui@genomics.cn

=head1 Version
    v1.0;    2012-5-21
    v1.1;    2012-6-6
    v1.3;    2013-5-22

=head1 NOTE
    if *split_fq.sh failed, clean $outdir and rerun this program.
=cut

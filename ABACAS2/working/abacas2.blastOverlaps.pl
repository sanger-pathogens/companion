#!/usr/bin/perl

# mh12@sanger.ac.uk

use strict;
use warnings;
use List::Util qw(min max);
use Getopt::Long;
use mh12_utils;

my %options = (
    annotate      => 0,
    annotate_only => 0,
    endTolerance  => 20,           # cutoff value of distance between start/end of read and start/end of overlap
    excludeFile   => "none",
    genome_length => 0,
    maxHitNo      => 5000,
    minLength     => 100,          # minimum length to count as an overlap
    nograph       => 0,
    range         => "85,95,99",

    #overhang       => 0,
    read_length     => 0,
    report_overlaps => 0,
    split           => 0,
    useOnlyFile     => "none",
    verbose         => 0,
);

my $ops_ok = GetOptions(
    \%options,     'annotate', 'annotate_only',  'endTolerance=i', 'excludeFile=s', 'genome_length=i', 'maxHitNo=i',
    'minLength=i', 'nograph',  'range=s',        'read_length=i', 'report_overlaps', 'useOnlyFile=s',
    'verbose',
);

if ($#ARGV < 3 or !($ops_ok)) {
    print
      "usage option 1:\n\tblastOverlaps.pl <normal> <ref fasta> <query fasta> <blast file (-m 8 or 9 output)> <output prefix> [options]\n\n",
      "usage option 2:\n\tblastOverlaps.pl <blast|blat> <ref fasta> <query fasta> <output prefix> [options]\n\n",
      "usage option 3:\n\tblastOverlaps.pl <combine> <ref fasta> <ref fasta split prefix> <ref_split prefix> <out prefix> [options]\n\n",
      "options:\n",
      "\t-annotate\tuse this to make an annotated blast file\n",
      "\t-annotate_only\tUse this to only make an annotated blast file, no other output\n",
      "\t-endTolerance <max dist to ends to count a hit> default: 20\n",
      "\t-excludeFile <file of readnames to ignore | none> default: none\n",
      "\t-genome_length <length in bp of genome> default: not used (i.e. size unknown)\n",
      "\t\tIf used, error plot of size/coverage will be made\n",
      "\t-maxHitNo <N> the max number of hits per read in blast file should be N (default=5000)\n",
      "\t-minLength <minimum hit length to count overlap> default: 100\n",
      "\t-nograph\tuse this to not plot graphs (will save histogram data, i.e. use this for split jobs)\n",
      "\t-range <list of comma-separated ranges> default: 85,95,99\n",
      "\t-read_length <N> use this for solexa, when all reads are the same length N (default=0 for unknown)\n",
      "\t-report_overlaps\tUse this to report the overlaps for each read (maybe a big file!)\n",
      "\t-useOnlyFile <file of readnames to use | none> default: none\n",
      "note: it is not recommended to use both the -exclude and -useOnly options.\n",
      "Using them will lead to paradoxes and unpredictable results\n\n";
    exit;
}

$options{annotate} = 1 if $options{annotate_only};

my $number_of_reads;
my %overlap_length_hist;
my %overlap_identity_hist;
my %hit_pc_identity_bins;
my %readlengths;           # to store read name->read length
my @overlap_count_hist;    # nth element is # read with n overlaps

$options{type} = $ARGV[0];

if ($ARGV[0] eq "normal") {
    $options{ref_fasta}  = $ARGV[1];
    $options{qry_fasta}  = $ARGV[2];
    $options{blast_in}   = $ARGV[3];
    $options{out_prefix} = $ARGV[4];
}
elsif ($ARGV[0] eq "blat" or $ARGV[0] eq "blast") {
    $options{ref_fasta}  = $ARGV[1];
    $options{qry_fasta}  = $ARGV[2];
    $options{out_prefix} = $ARGV[3];
}
elsif ($ARGV[0] eq "combine") {
    $options{ref_fasta}              = $ARGV[1];
    $options{ref_fasta_split_prefix} = $ARGV[2];
    $options{split_prefix}           = $ARGV[3];
    $options{out_prefix}             = $ARGV[4];
}
else {
    die "Confusion error: unknown option $ARGV[0] used\n\n";
}

print "Getting the read lengths..." if $options{verbose};

# get the read lengths
$number_of_reads = fill_readlengths_hash(\%options, \%readlengths);

#$number_of_reads = 5320094;
#$readlengths{all_reads_length} = 76;

print "$number_of_reads reads found\n" if $options{verbose};

if ($options{type} ne "combine") {
    findOverlaps(\%options, \%overlap_length_hist, \%overlap_identity_hist, \%hit_pc_identity_bins, \%readlengths,
        \@overlap_count_hist, $number_of_reads);
}
else {
    print "checking split files to be combined..." if $options{verbose};

    # check if it looks like all the split jobs ran ok
    my @ls = `ls -1 $options{split_prefix}.[0-9]* $options{ref_fasta_split_prefix}.[0-9]*`;
    chomp @ls;
    my $count1 = 0;
    my $count2 = 0;

    foreach (@ls) {
        if    (/$options{ref_fasta_split_prefix}\.\d+$/)            {$count1++}
        elsif (/$options{split_prefix}\.\d+.*overlap_length.*gram/) {$count2++}
    }

    if ($count1 == $count2 || $count1 != 0) {
        print "found $count1 split jobs to be combined\n";
    }
    else {
        die "Looks like something is wrong with the split jobs...expected:$count1 found:$count2\n";
    }

    # combine the hits_per_read_info files
    system
      "zcat $options{split_prefix}.[0-9]*.hits_per_read_info.gz | gzip -9 > $options{out_prefix}.hits_per_read_info.gz";

    print "hits_per_read_info files combined\nCombining overlap_length files..." if $options{verbose};

    # combine the overlap_length.histogram files into a hash
    my @data =
      `paste $options{split_prefix}.[0-9]*.overlap_length.histogram | awk '{s=0; for (i=2; i<=NF; i+=2) {s+=\$i}; print \$1,s }'`;
    chomp @data;
    array2hash(\@data, \%overlap_length_hist);

    print "done\nCombining overlap_identity files..." if $options{verbose};

    # combine the overlap_identity.histogram files into a hash
    @data =
      `paste $options{split_prefix}.[0-9]*.overlap_identity.histogram | awk '{s=0; for (i=2; i<=NF; i+=2) {s+=\$i}; print \$1,s }'`;
    chomp @data;
    array2hash(\@data, \%overlap_identity_hist);

    print "done\n" if $options{verbose};

    # if present, combine the annotated bla(s)t files
    if (-e "$options{split_prefix}.1.blast.annotated.gz") {
        system
          "zcat $options{split_prefix}.[0-9]*.blast.annotated.gz | gzip -9 > $options{out_prefix}.blast.annotated.gz";
        print "Annotated blast files found and combined\n" if $options{verbose};
    }

    # if present, combine the overlap report files
    if (-e "$options{split_prefix}.1.overlap_report.gz") {
        system
          "zcat $options{split_prefix}.[0-9]*.overlap_report.gz | gzip -9 > $options{out_prefix}.overlap_report.gz";
        print "Overlap report files found and combined\n" if $options{verbose};
    }

    # combine the overlap_count.histogram files into an array
    print "Combining overlap_count files..." if $options{verbose};
    @data = `paste $options{split_prefix}.[0-9]*.overlap_count.histogram`;
    chomp @data;

    # find out how many ranges were given
    my @data_line = split /\s+/, $data[0];
    my $range_count = 0;

    if ($data_line[0] eq "read") {
        foreach (@data_line[1 .. $#data_line]) {
            last if ($_ eq "read");
            $range_count++;
        }
    }
    else {
        die "format of overlap_identity.histogram files not recognised\nThis shouldn't have happened\n\n";
    }

    $options{range} = join ",", @data_line[1 .. $range_count];

    print "ranges found: $options{range}..." if $options{verbose};

    # combine the results into the overlap_count array
    foreach (@data[1 .. $#data]) {
        @data_line = split /\W+/;
        my @a;
        foreach (1 .. $range_count) {push @a, 0;}

        foreach my $i (1 .. $#data_line) {
            $i % 4 == 0 or $a[($i % 4) - 1] += $data_line[$i];
        }

        push @overlap_count_hist, [@a];

    }

    print "done\n" if $options{verbose};
}

exit if $options{annotate_only};

print "Writing histograms and plot files..." if $options{verbose};

write_hists_and_graphs(\%options, \%overlap_length_hist, \%overlap_identity_hist, \%hit_pc_identity_bins, \%readlengths,
    \@overlap_count_hist, $number_of_reads);

print "done\n" if $options{verbose};

#______________________________________________________________________________
#
#                            SUBROUTINES
#______________________________________________________________________________
sub array2hash {
    my $a = shift;
    my $h = shift;

    foreach (@$a) {
        my @tmp = split /\s/;
        $h->{ $tmp[0] } = $tmp[1];
    }
}


#usage: fill_readlengths_hash(ops_hash_ref, realength_hash_ref)
#returns: number of reads
sub fill_readlengths_hash {
    my $ops     = shift;
    my $lengths = shift;
    my $no_of_reads;

    # fill the readlength hash if readlength not known
    if ($ops->{read_length}) {
        $lengths->{all_reads_length} = $ops->{read_length};
        $no_of_reads = `grep -c ">" $ops->{ref_fasta}`;
        chomp $no_of_reads;
    }
    else {
        mh12_utils::fasta2length($ops->{ref_fasta}, $lengths);

        if ($ops->{qry_fasta}){
            mh12_utils::fasta2length($ops->{qry_fasta}, $lengths);
        }

        $no_of_reads = scalar keys %$lengths;
    }

    return $no_of_reads;
}

sub findOverlaps {
    my $ops_ref                   = shift;
    my $overlap_length_hist_ref   = shift;
    my $overlap_identity_hist_ref = shift;
    my $hit_pc_identity_bins_ref  = shift;
    my $readlengths_ref           = shift;
    my $overlap_count_hist_ref    = shift;
    my $no_of_reads               = shift;

    my @ranges = split /,/, $ops_ref->{range};

    my $overlap_length;    # to store the length of each overlap > minLength
    my $current_readname = "none";
    my @current_percents;
    my %exclude;           # stores the reads to exclude
    my %useOnly;           # stores the reads to only use
    my $overlap_report = "";
    my $blast_filehandle;
    my $annotate_outstring = "";                # stores the annotated blast file output
    my @array_of_zeros     = (0 .. $#ranges);

    foreach (@array_of_zeros) {
        $array_of_zeros[$_] = 0;
    }

    foreach (0 .. $ops_ref->{maxHitNo}) {
        push @{$overlap_count_hist_ref}, [@array_of_zeros];
    }

    # fill the excluded reads/ useOnly hash(es)
    if ($ops_ref->{excludeFile} ne "none") {file2hash($ops_ref->{excludeFile}, \%exclude)}
    if ($ops_ref->{useOnlyFile} ne "none") {file2hash($ops_ref->{useOnlyFile}, \%useOnly)}

    # find the longest read length
    my $longest_read_length = 0;

    if ($ops_ref->{read_length}) {
        $longest_read_length = $ops_ref->{read_length};
    }
    else {
        $_ > $longest_read_length and $longest_read_length = $_ for values %$readlengths_ref;
    }

    # initialise hashes which will store the histograms of overlap lengths
    # and overlap identities
    for (-$longest_read_length .. $longest_read_length) {
        $overlap_length_hist_ref->{$_} = 0;
    }

    #for (my $i=$ranges[0]; $i <= 100; $i += 0.1){
    for (-1000 .. 1000) {

        #$i = sprintf("%.2f", $i);    # round $i to the nearest tenth
        #$overlap_identity_hist_ref->{ sprintf("%.2f", $_) } = 0;    #rounded to nearest 1/10th
        $overlap_identity_hist_ref->{ sprintf("%.1f", $_ / 10) } = 0;    #rounded to nearest 1/10th
    }

    print "Getting query readnames..." if $ops_ref->{verbose};

    # initialise hash (some reads might not be in the blast file)
    my @readnames = `awk '/^>/ {print substr(\$1,2,length(\$1)-1)}' $options{qry_fasta}`;
    my @zeros     = ();

    for my $i (0 .. $#ranges) {
        push @zeros, 0;
    }

    foreach (@readnames) {
        chomp;
        $hit_pc_identity_bins_ref->{$_} = [@zeros];
    }

    undef @readnames;

    print "Number of reads found: " . (scalar keys %{$hit_pc_identity_bins_ref}) . "\n" if $ops_ref->{verbose};

    if ($options{type} eq "normal") {

        # check if blast file is gzipped
        if ($ops_ref->{blast_in} =~ /\.gz$/) {
            open $blast_filehandle, "zcat $ops_ref->{blast_in} && echo \"lastline\" |"
              or die "error opening zipped file $ops_ref->{blast_in}\n";
        }
        else {
            open $blast_filehandle, "cat $ops_ref->{blast_in} && echo \"lastline\" |"
              or die("Error opening $ops_ref->{blast_in}");
        }
    }
    elsif ($options{type} eq "blat") {
        print "Running blat..." if $ops_ref->{verbose};

        open $blast_filehandle,
          "blat $ops_ref->{ref_fasta} $ops_ref->{qry_fasta} stdout -out=blast8 && echo \"lastline\" |"
          or die "error running blat on $ops_ref->{ref_fasta}\n";

        print "done\n" if $ops_ref->{verbose};
    }
    elsif ($options{type} eq "blast") {
        print "Running blast..." if $ops_ref->{verbose};
        unless (-e "$ops_ref->{ref_fasta}.nhr") {
            system "formatdb -p F -i $ops_ref->{ref_fasta}";
        }

        open $blast_filehandle,
          "blastall -p blastn -W 15 -F F -e 1e-10  -d $ops_ref->{ref_fasta} -i $ops_ref->{qry_fasta} -m 8 && echo \"lastline\" | "
          or die "error running blast on $ops_ref->{ref_fasta}\n";

        print "done\n" if $ops_ref->{verbose};
    }

    print "Parsing bla(s)t results..." if $ops_ref->{verbose};

    # main loop of file - go through each line of blast file and decide whether
    # each one is an overlap or not
    while (my $line = <$blast_filehandle>) {

        # check if line is a blast hit or a ssaha hit and put $1=query name, $2=subject name
        if ($line =~ /^ALIGNMENT::\d+\s+\S+\s+(\S+)\s+(\S+)/ || $line =~ /(^[^\#]\S+)\s(\S+)/ || $line =~ /^lastline$/)
        {
            chomp $line;

            # skip this read?
            if (
                !($line eq "lastline")
                && (   $exclude{$1}
                    || $exclude{$2}
                    || ($ops_ref->{useOnlyFile} ne "none" && !($useOnly{$1} && $useOnly{$2})))
              )
            {

                # print line (if wanted)
                if ($ops_ref->{annotate}) {
                    $annotate_outstring .= "$line EXCLUDED\n";
                }

                next;
            }

            # if 1st time round, set current_readname
            if ($current_readname eq "none") {
                $current_readname = $1;
            }
            else {
                if ($line eq "lastline" || $current_readname ne $1) {    # if current name != next name
                                                                         # reset the counts array
                    my @count = ();
                    for my $i (0 .. (1 + 2 * $#ranges)) {$count[$i] = 0;}

                    # count up the hits in each % range - put each hit in a bin depending
                    # on its value.  e.g. for range values r1, r2, r3, there will be six
                    # bins: three for overlaps (pc>0) and three for overhangs (pc<0):
                    foreach my $pc (@current_percents) {
                        for my $i (0 .. $#ranges) {
                            if ($pc > 0 && $pc >= $ranges[$#ranges - $i]) {
                                $count[$#ranges - $i]++;
                                last;
                            }
                            elsif ($pc < 0 && -$pc >= $ranges[$#ranges - $i]) {
                                $count[2 * $#ranges - $i]++;
                                last;
                            }
                        }
                    }

                    # there may be more hits than the maximum - check for this and chop
                    # them down to maxHitNo
                    foreach (@count) {
                        if ($_ > $ops_ref->{maxHitNo}) {
                            $_ = $ops_ref->{maxHitNo};
                        }
                    }

                    $hit_pc_identity_bins_ref->{$current_readname} = [@count];

                    # print "@count $count[$#ranges]\n";
                    foreach (0 .. $#ranges) {

                        #$overlap_count_hist_ref->[$_][$count[$_]]++;
                        $overlap_count_hist_ref->[$count[$_]][$_]++;
                    }

                    #$overlap_count_hist_ref->[$count[$#ranges]]++;

                    if ($line eq "lastline") {last}

                    $current_readname = $1;    # reset current readname
                    @current_percents = ();    # reset percentage array
                }
            }

            # get the %idenity, if the read is an overlap (=0 otherwise)
            my $percent = isOverlap(\$line, $ops_ref->{minLength}, $ops_ref->{endTolerance},
                $readlengths_ref, $overlap_length_hist_ref);

            # add %identity to ouptut
            if ($percent != 0) {
                push(@current_percents, $percent);

                #my $p = sprintf("%.2f", $percent);    # round $i to the nearest tenth
                #$overlap_identity_hist_ref->{ sprintf("%.2f", $percent) }++;
                $overlap_identity_hist_ref->{ sprintf("%.1f", $percent) }++;
            }



            my %h;
            my $lengths;

            if ( $ops_ref->{annotate} || ($ops_ref->{report_overlaps} && $percent > 0)){
                if ($ops_ref->{read_length}){
                    $lengths = "$ops_ref->{read_length} $ops_ref->{read_length}";
                }
                else{
                   blast_or_ssaha2overlap_hash(\$line, \%h);
                   $lengths = "$readlengths{$h{q_name}} $readlengths{$h{s_name}}";
                }
            }

            # print annotated blast line if wanted
            if ($ops_ref->{annotate}) {
                # add either "OVERLAP", "OVERHANG", or "NO OVERLAP" to the annotated blast file line
                if ($percent > 0) {
                    $annotate_outstring .= "$line $lengths OVERLAP\n";
                }
                elsif ($percent < 0) {
                    $annotate_outstring .= "$line $lengths OVERHANG\n";
                }
                else {
                    $annotate_outstring .= "$line $lengths NA\n";
                }
            }

            if ($ops_ref->{report_overlaps} && $percent > 0) {
                $overlap_report .= "$line $lengths\n";
            }
        }
    }

    close F;

    print "done\n" if $ops_ref->{verbose};

    if ($ops_ref->{report_overlaps}) {
        open F, "| gzip -9 > $ops_ref->{out_prefix}.overlap_report.gz"
          or die "error opening $ops_ref->{out_prefix}.overlap_report.gz\n";
        print F $overlap_report;
        close F;
    }

    # write the annotated blast file
    if ($ops_ref->{annotate}) {
        open F, "| gzip -9 > $ops_ref->{out_prefix}.blast.annotated.gz"
          or die "Error opening $ops_ref->{out_prefix}.blast.annotated\n";
        print F $annotate_outstring;
        close F;
    }
}

#______________________________________________________________________________
sub write_hists_and_graphs {
    my $ops_ref                   = shift;
    my $overlap_length_hist_ref   = shift;
    my $overlap_identity_hist_ref = shift;
    my $hit_pc_identity_bins_ref  = shift;
    my $readlengths_ref           = shift;
    my $overlap_count_hist_ref    = shift;
    my $no_of_reads               = shift;

    my @ranges      = split /,/, $ops_ref->{range};
    my @coverage    = (0 .. $ops_ref->{maxHitNo});    # coverage[n] = predicted coverage with cutoff=n
    my @genome_size = (0 .. $ops_ref->{maxHitNo});    # genome_size[n] = predicted genome size with cutoff=n
    my $mean_read_length;
    my @R_overlap_lengths;                            # stores overlap lengths for R plotting later

    # write the file of number-of-overlaps-each-read-has histogram
    open F, ">$ops_ref->{out_prefix}.overlap_count.histogram"
      or die "error opening $ops_ref->{out_prefix}.overlap_count.histogram\n";
    print F "read @ranges\n";

    foreach my $i (0 .. $ops_ref->{maxHitNo}) {
        print F "$i ";
        foreach my $j (0 .. $#ranges) {
            print F "$overlap_count_hist_ref->[$i][$j] ";
        }
        print F "\n";
    }
    close F;

    # write the file of overlap lengths histogram
    open F, ">$ops_ref->{out_prefix}.overlap_length.histogram"
      or die "error opening $ops_ref->{out_prefix}.overlap_length.histogram\n";

    foreach (sort {$a <=> $b} keys %$overlap_length_hist_ref) {
        print F "$_ $overlap_length_hist_ref->{$_}\n";

        if ($_ >= 0) {
            push @R_overlap_lengths, $_;
        }
    }
    close F;

    # write the file of %identity histogram
    open F, ">$ops_ref->{out_prefix}.overlap_identity.histogram"
      or die "error opening $ops_ref->{out_prefix}.overlap_identity.histogram\n";

    foreach (sort {$a <=> $b} keys %$overlap_identity_hist_ref) {
        print F "$_ $overlap_identity_hist_ref->{$_}\n";
    }
    close F;

    # get the mean read length
    if ($ops_ref->{read_length}) {    # if all reads the same length
        $mean_read_length = $ops_ref->{read_length};
    }
    else {                            # otherwise, need to calculate the mean read length
        my $sum = 0;

        foreach (keys %{$readlengths_ref}) {
            $sum += $readlengths_ref->{$_};
        }

        $mean_read_length = $sum / (scalar keys %{$readlengths_ref});

    }

    open F, ">$ops_ref->{out_prefix}.genome_size_etc.plot"
      or die "error opening $ops_ref->{out_prefix}.genome_size_etc.plot\n";

    # work out the genome size and coverage for each cutoff value
    my @mean_overlap = (0 .. $ops_ref->{maxHitNo});

    #my $cumulative_count    = $overlap_count_hist_ref->[$#ranges][0];
    my $cumulative_count    = $overlap_count_hist_ref->[0][$#ranges];
    my $cumulative_overlaps = 0;

    foreach my $i (1 .. $ops_ref->{maxHitNo}) {

        #$cumulative_count += $overlap_count_hist_ref->[$#ranges][$i];
        #$cumulative_overlaps += $i * $overlap_count_hist_ref->[$#ranges][$i];
        $cumulative_count += $overlap_count_hist_ref->[$i][$#ranges];
        $cumulative_overlaps += $i * $overlap_count_hist_ref->[$i][$#ranges];

        if ($cumulative_count > 0){
            $mean_overlap[$i] = 1.0 * $cumulative_overlaps / $cumulative_count;
        }
        else {
            $mean_overlap[$i] = 0;
        }

        $coverage[$i] = 0.5 * $mean_overlap[$i] * $mean_read_length / ($mean_read_length - $ops_ref->{minLength});

        if ($coverage[$i] > 0){
            $genome_size[$i] = $cumulative_count * $mean_read_length / $coverage[$i];
        }
        else {
            $genome_size[$i] = 0
        }


        print F "$i $coverage[$i] $genome_size[$i]";

        if ($ops_ref->{genome_length}) {
            print F " " . (100 * ($ops_ref->{genome_length} - $genome_size[$i]) / $ops_ref->{genome_length});
        }

        print F "\n";

        #print "cutoff: $i, mean:$mean_overlap[$i], c:$coverage[$i], G:$genome_size[$i]\n";
    }

    close F;
    unless ($ops_ref->{type} eq "combine") {

        # write the file of number of hits each read has
        open F, "| gzip -9 >$ops_ref->{out_prefix}.hits_per_read_info.gz"
          or die "Error opening $ops_ref->{out_prefix}.hits_per_read.info.gz\n";

        foreach (sort keys %$hit_pc_identity_bins_ref) {
            print F "$_ " . join(" ", @{ $hit_pc_identity_bins_ref->{$_} }) . "\n";
        }

        close F;
    }

    # make an R script to make the plots
    unless ($ops_ref->{nograph}) {
        open F, ">$ops_ref->{out_prefix}.make_plots.R" or die "error opening $ops_ref->{out_prefix}.make_plots.R\n";

        #________________________ make overlap length histogram __________________________
        print F 'x=c(';
        print F join ',', @R_overlap_lengths;
        print F ")\n";

        print F 'y=c(';
        my $tmpstring = "";
        foreach (@R_overlap_lengths) {
            $tmpstring .= "$overlap_length_hist_ref->{$_},";
        }
        print F substr($tmpstring, 0, -1);
        print F ")\n",
          "pdf(\"$ops_ref->{out_prefix}.overlap_length.histogram.pdf\")\n",
          "    plot(x,y,xlab=\"Overlap Length\",ylab=\"Frequency\",type=\"l\")\n",
          "dev.off()\n\n\n";

        #_____________________________ make overlap identity histogram ___________________
        $tmpstring = "x=c(";
        my $tmpstring2 = "y=c(";
        my $flag       = 0;

        foreach (sort {$a <=> $b} keys %$overlap_identity_hist_ref) {
            if ($_ >= 85) {
                if ($flag || $overlap_identity_hist_ref->{$_} > 0) {
                    $tmpstring  .= "$_,";
                    $tmpstring2 .= "$overlap_identity_hist_ref->{$_},";
                    $flag = 1;
                }
            }
        }

        print F substr($tmpstring,  0, -1) . ")\n";
        print F substr($tmpstring2, 0, -1) . ")\n";

        # make labels for y axis: absolute count, %
        print F "percent  =  round( c(1:5) * max(y) / 5 )\n",
          "tick_pos = signif( c(1:5) * max(y) / 5, 2 )\n",
          "ylabels = paste(percent, \"\%\", sep=\"\")\n",
          "pdf(\"$ops_ref->{out_prefix}.overlap_identity.histogram.pdf\")\n",
          "   par(mar=c(5, 4, 4, 5) + 0.1)\n",
          "   barplot(y,x, xlab=\"Overlap Percent Identity\", ylab=\"Frequency\", names.arg=x)\n",
          "   axis(4, at=c(0, tick_pos), labels=c(0, ylabels))\n",
          "   mtext(\"Percentage of Reads\", side=4, line=3, cex.lab=1,las=3)\n",
          "dev.off()\n\n\n";

        #___________________ make plot of genome size/coverage vs cutoff __________________
        print F "cov=c(" . (join ",", @coverage[1 .. $#coverage]) . ")\n",
          "gsize=c(" . (join ",", @genome_size[1 .. $#genome_size]) . ")\n",
          "x=c(1:$#coverage)\n";

        # make tick mark names vectors
        print F "maxc = max(cov); minc = min(cov)\n", "maxg = max(gsize); ming = min (gsize)\n";

        # if coverage and genome size remain constant, then min=max, which will
        # cause a divide by zero error later.  Hence scale them a little if this happens
        print F "if(minc == maxc){minc = minc * 0.9; maxc = maxc * 1.1}\n",
          "if(ming == maxg){ming = ming * 0.9; maxg = maxg * 1.1}\n";

        print F "tick_pos = 0.25 * c(0:4)\n",
          "tick_labelsc = signif(0.25 * c(0:4) * (maxc-minc) + minc, 2)\n",
          "tick_labelsg =  round(0.25 * c(0:4) * (maxg-ming) + ming)\n";

        # scale the coverage/genome size values in d to the range [0,1]
        print F "cov = ( cov - minc ) / ( maxc - minc )\n", "gsize = ( gsize - ming ) / ( maxg - ming )\n";

        print F "pdf(\"$ops_ref->{out_prefix}.genome_size_coverage_plot.pdf\")\n",
          "    par(mar=c(4,4,5,5))\n",
          "    plot (x,cov, type=\"l\", col=\"red\", xlab=\"Cutoff\", yaxt=\"n\", ylab=\"\")\n",
          "    lines(x,gsize, col=\"blue\")\n",
          "    axis(2, at=c(tick_pos), labels=c(tick_labelsc))\n",
          "    axis(4, at=c(tick_pos), labels=c(tick_labelsg))\n",
          "    mtext(\"Coverage\", side=2, line=3, cex.lab=1, col=\"red\")\n",
          "    mtext(\"Genome Size\", side=4, line=3, cex.lab=1, col=\"blue\")\n",
          "dev.off()\n\n\n";

        #___________________ make plot of genome size/coverage vs cutoff > 10 __________________
        # do the same plot, but with the first 10 values cutoff (they can be massively different
        # which makes the plot look a bit rubbish)
        print F "cov=c(" . (join ",", @coverage[10 .. $#coverage]) . ")\n",
          "gsize=c(" . (join ",", @genome_size[10 .. $#genome_size]) . ")\n",
          "x=c(10:$#coverage)\n";

        # make tick mark names vectors
        print F "maxc = max(cov); minc = min(cov)\n", "maxg = max(gsize); ming = min (gsize)\n";

        # if coverage and genome size remain constant, then min=max, which will
        # cause a divide by zero error later.  Hence scale them a little if this happens
        print F "if(minc == maxc){minc = minc * 0.9; maxc = maxc * 1.1}\n",
          "if(ming == maxg){ming = ming * 0.9; maxg = maxg * 1.1}\n";

        print F "tick_pos = 0.25 * c(0:4)\n",
          "tick_labelsc = signif(0.25 * c(0:4) * (maxc-minc) + minc, 2)\n",
          "tick_labelsg =  round(0.25 * c(0:4) * (maxg-ming) + ming)\n";

        # scale the coverage/genome size values in d to the range [0,1]
        print F "cov = ( cov - minc ) / ( maxc - minc )\n", "gsize = ( gsize - ming ) / ( maxg - ming )\n";

        print F "pdf(\"$ops_ref->{out_prefix}.genome_size_coverage_plot-over10.pdf\")\n",
          "    par(mar=c(4,4,5,5))\n",
          "    plot (x,cov, type=\"l\", col=\"red\", xlab=\"Cutoff\", yaxt=\"n\", ylab=\"\")\n",
          "    lines(x,gsize, col=\"blue\")\n",
          "    axis(2, at=c(tick_pos), labels=c(tick_labelsc))\n",
          "    axis(4, at=c(tick_pos), labels=c(tick_labelsg))\n",
          "    mtext(\"Coverage\", side=2, line=3, cex.lab=1, col=\"red\")\n",
          "    mtext(\"Genome Size\", side=4, line=3, cex.lab=1, col=\"blue\")\n",
          "dev.off()\n\n\n";




        #________________________ if genome size known, make error plot __________________
        if ($ops_ref->{genome_length}) {
            print F "x=c(1:$#coverage)\n";
            $tmpstring = "y=c(";

            foreach (@genome_size[1 .. $#genome_size]) {
                $tmpstring .= (100 * ($_ - $ops_ref->{genome_length}) / $ops_ref->{genome_length}) . ",";
            }

            print F substr($tmpstring, 0, -1) . ")\n",
              "pdf(\"$ops_ref->{out_prefix}.genome_size_error.pdf\")\n",
              "    plot (x,y, type=\"l\", xlab=\"Cutoff\", ylab=\"Genome size \% error\")\n",
              "dev.off()\n\n\n";

            $tmpstring = "y=c(";


            print F "x=c(10:$#coverage)\n";
            foreach (@genome_size[10 .. $#genome_size]) {
                $tmpstring .= (100 * ($_ - $ops_ref->{genome_length}) / $ops_ref->{genome_length}) . ",";
            }

            print F substr($tmpstring, 0, -1) . ")\n",
              "pdf(\"$ops_ref->{out_prefix}.genome_size_error-over10.pdf\")\n",
              "    plot (x,y, type=\"l\", xlab=\"Cutoff\", ylab=\"Genome size \% error\")\n",
              "dev.off()\n\n\n";
        }

        #_______________________ make histogram of number of overlaps each read has__________________
        # load in the histogram data
        print F "d = read.table(\"$ops_ref->{out_prefix}.overlap_count.histogram\", header=TRUE)\n";

        # make tick marks and labels for % of total reads for y axis
        print F "tick_pos  = c(1:5) * max(d[,2:dim(d)[2]]) / 5\n",
          "percent =    round( 100 * tick_pos / $no_of_reads )\n";

        # construct the plot legend
        print F "plot_legend = vector(length=dim(d)[2]-1)\n",
          "l = length(plot_legend)\n",
          "for (i in 1:(l-1)){\n",
          "    plot_legend[i] = paste(substr(colnames(d)[i+1], 2, 5), \"\%-\", substr(colnames(d)[i+2], 2, 5), \"\% identity\", sep=\"\")\n",
          "}\n",
          "plot_legend[l] = paste(\">\", substr(colnames(d)[l+1], 2, 5),  \"\% identity\", sep=\"\")\n";

        # make the barplot
        print F "pdf(\"$ops_ref->{out_prefix}.overlap_count.histogram.pdf\")\n",
          "  par(mar=c(4,4,5,5))\n",
          "  barplot(t(as.matrix(d[1:50,2:dim(d)[2]])), col=rainbow(dim(d)[2] - 1), legend=plot_legend, beside=TRUE, names.arg=d[1:50,1], xlab=\"Number of overlaps\", ylab=\"Frequency\")\n",
          "  axis(4, at=c(0,tick_pos), labels=c(0,percent))\n",
          "  mtext(\"Percent of total reads\", side=4, line=3)\n",
          "dev.off()\n";

        # make plots of >99% id hits of larger values on x axis to see repeats
        print F "d10 = d[10:dim(d)[1], dim(d)[2]]\n",
          "d20 = d[20:dim(d)[1], dim(d)[2]]\n",
          "d30 = d[30:dim(d)[1], dim(d)[2]]\n",
          "d40 = d[40:dim(d)[1], dim(d)[2]]\n",
          "pdf(\"$ops_ref->{out_prefix}.overlap_count.histogram-over10.pdf\", width=12, height=6)\n",
          "    par(mfrow=c(2,2))\n",
          "    barplot(d10, names.arg=d[10:dim(d)[1],1], xlab=\"Number of overlaps\", ylab=\"Frequency\")\n",
          "    barplot(d20, names.arg=d[20:dim(d)[1],1], xlab=\"Number of overlaps\", ylab=\"Frequency\")\n",
          "    barplot(d30, names.arg=d[30:dim(d)[1],1], xlab=\"Number of overlaps\", ylab=\"Frequency\")\n",
          "    barplot(d40, names.arg=d[40:dim(d)[1],1], xlab=\"Number of overlaps\", ylab=\"Frequency\")\n",
          "dev.off()\n";

        close F;

        # run the R script
        system "R CMD BATCH $ops_ref->{out_prefix}.make_plots.R";
    }
}

sub blast_or_ssaha2overlap_hash {
    my $in_ref     = shift;
    my $out_ref = shift;

    my @tmp = split /\s+/, $$in_ref;

    if ($$in_ref =~ /ALIGNMENT/) {    # if the inputline is a ssaha alignment
        $out_ref->{q_name}  = $tmp[2];       # query nameA
        $out_ref->{s_name}  = $tmp[3];       # subject name
        $out_ref->{q_start} = $tmp[4];       # query start
        $out_ref->{q_end}   = $tmp[5];       # query end
        $out_ref->{s_start} = $tmp[6];       # subject start
        $out_ref->{s_end}   = $tmp[7];       # subject end
        $out_ref->{pc_id}   = $tmp[10];      # % identity
    }
    else {                              # else, it's a blast line
        $out_ref->{q_name}  = $tmp[0];       # query nameA
        $out_ref->{s_name}  = $tmp[1];       # subject name
        $out_ref->{q_start} = $tmp[6];       # query start
        $out_ref->{q_end}   = $tmp[7];       # query end
        $out_ref->{s_start} = $tmp[8];       # subject start
        $out_ref->{s_end}   = $tmp[9];       # subject end
        $out_ref->{pc_id}   = $tmp[2];       # % identity
    }
}



#______________________________________________________________________________
# takes line of output from blast, returns the %identity if it is an overlap of
# two reads, -1 * percent identity if it is an overang, else returns 0.
sub isOverlap {
    my $inputline_ref       = shift;
    my $length_cutoff   = shift;
    my $end_range       = shift;
    my $length_hash_ref = shift;
    my $overlap_ref     = shift;

    my %data;    # will store the relevant data from the input line

    blast_or_ssaha2overlap_hash($inputline_ref, \%data);

#    if ($inputline =~ /ALIGNMENT/) {    # if the inputline is a ssaha alignment
#        $data{q_name}  = $tmp[2];       # query nameA
#        $data{s_name}  = $tmp[3];       # subject name
#        $data{q_start} = $tmp[4];       # query start
#        $data{q_end}   = $tmp[5];       # query end
#        $data{s_start} = $tmp[6];       # subject start
#        $data{s_end}   = $tmp[7];       # subject end
#        $data{pc_id}   = $tmp[10];      # % identity
#    }
#    else {                              # else, it's a blast line
#        $data{q_name}  = $tmp[0];       # query nameA
#        $data{s_name}  = $tmp[1];       # subject name
#        $data{q_start} = $tmp[6];       # query start
#        $data{q_end}   = $tmp[7];       # query end
#        $data{s_start} = $tmp[8];       # subject start
#        $data{s_end}   = $tmp[9];       # subject end
#        $data{pc_id}   = $tmp[2];       # % identity
#    }

    # if the hit is a read to itself,then don't count it as an overlap
    if ($data{q_name} eq $data{s_name}) {return 0;}

    # get the start/end values of hit
    my $s1 = $data{q_start};
    my $e1 = $data{q_end};
    my $s2 = $data{s_start};
    my $e2 = $data{s_end};

    my $q_length;
    my $s_length;

    if ($length_hash_ref->{all_reads_length}) {
        $q_length = $length_hash_ref->{all_reads_length};
        $s_length = $length_hash_ref->{all_reads_length};
    }
    else {
        $q_length = $length_hash_ref->{ $data{q_name} };
        $s_length = $length_hash_ref->{ $data{s_name} };
    }

    # reorder the starts/ends so that start <end
    if ($s1 > $e1) {
        $s1 = $q_length - $data{q_start} + 1;
        $e1 = $q_length - $data{q_end} + 1;
    }

    if ($s2 > $e2) {
        $s2 = $s_length - $data{s_start} + 1;
        $e2 = $s_length - $data{s_end} + 1;
    }

    # if the start/end of an overlap is within $end_range of the end of the sequence,
    # then assume the overlap starts at the start/end
    if ($s1 < $end_range)             {$s1 = 1}
    if ($s2 < $end_range)             {$s2 = 1}
    if ($e1 > $q_length - $end_range) {$e1 = $q_length}
    if ($e2 > $s_length - $end_range) {$e2 = $s_length}

    # if the length of the overlap is < cutoff, don't count it
    if ($e1 - $s1 + 1 < $length_cutoff) {return 0;}

    my $test_s1 = ($s1 == 1);
    my $test_e1 = ($e1 == $q_length);
    my $test_s2 = ($s2 == 1);
    my $test_e2 = ($e2 == $s_length);

    # if it is an overlap, add overlap length to legnth array, return the %identity
    if (($test_s1 && $test_e1) || ($test_s1 && $test_e2) || ($test_s2 && $test_e1) || ($test_s2 && $test_e2)) {
        my $l = ($e1 - $s1 + 1);
        $overlap_ref->{ ($e1 - $s1 + 1) }++;
        return $data{pc_id};
    }
    elsif (($test_s1 && $test_s2) || ($test_e1 && $test_e2)) {    # if it's an overhang
        $overlap_ref->{ -($e1 - $s1 + 1) }++;
        return -$data{pc_id};
    }
    else {                                                        # otherwise, not an overlap, so return 0
        return 0;
    }
}

#______________________________________________________________________________

# takes a file and fills a hash with: "line" => 1, for each line of file
#
# arg0: filename
# arg1: ref to hash
sub file2hash {
    my $file     = shift;
    my $hash_ref = shift;

    open F, $file or die("Error opening $file");
    while (my $line = <F>) {
        chomp $line;
        $hash_ref->{$line} = 1;
    }
}

#______________________________________________________________________________

# usage: array2histogram(array ref, bar width, outfile name);
sub array2histogram {
    my $array_ref = shift;
    my $width     = shift;
    my $outfile   = shift;
    my $scale     = shift;

    my $min  = min(@$array_ref);
    my $max  = max(@$array_ref);
    my %bins = ();

    # round min and max down to nearest tenth
    $min = $width * int($min / $width);
    $max = $width * int($max / $width);

    # fill hash with keys => 0
    for (my $i = $min; $i <= $max; $i += $width) {
        $i = sprintf("%.2f", $i);    # round $i to the nearest tenth
        $bins{$i} = 0;
    }

    my @keys = reverse sort {$a <=> $b} keys %bins;

    # count up the % identity values in each bin
    foreach my $val (@$array_ref) {
        foreach (@keys) {
            if ($val >= $_) {
                $bins{$_}++;
                last;
            }
        }
    }

    # scale the values, if necessary
    if ($scale != 0) {
        foreach (@keys) {
            $bins{$_} = 100 * $bins{$_} / $scale;

        }
    }

    # write the percent identity histogram file
    open F, ">$outfile" or die "Error opening $outfile";
    foreach (sort {$a <=> $b} keys %bins) {
        print F "$_ $bins{$_}\n";
        next;
    }
    close F;
}

#______________________________________________________________________________

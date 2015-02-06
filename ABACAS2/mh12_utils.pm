package mh12_utils;


use strict;
use warnings;
use Carp;

use base 'Exporter';

#our @EXPORT = qw(fasta2hash
#                 fasta2length
#                 file2hash
#                 read2mate
#                 read2plate);
#

# takes a fasta file and hash ref.  Fills the hash with id=>seq
# usage: fasta2length(fasta filename, hash reference)
sub fasta2hash{
    my $infile = shift;
    my $h_ref  = shift;

    my $id;
    my $seq;

    open F, $infile or croak "error opening $infile\n";
    my $s = join "", <F>;
    $s= substr($s,1);         # remove the first >
    $s =~ s/\n{2,}/\n/g;      # remove empty lines
    my @a = split /\n\>/, $s; # each element is id\nseq, where seq may have \n's
    undef $s;                 # free up some memory
    close F;

    # fill the hash
    foreach (@a){
        if (/(.*?\n)(.*)/s){
            $id = $1;
            chomp $id;

            $seq = $2;
            $seq =~ s/\n//g;

            $h_ref->{$id} = $seq;
        }
    }
}
#______________________________________________________________________________


# Takes a fasta file and hash reference. Fills the hash with id=>(array with
# 0 element: bases, 1st element: qual scores.  Fakes the quality scores to
# be all the same
#
# usage: fasta2fastq_hash(fasta filename, hash reference, [phred score (default:50)])
# (i.e. third parameter is optional: if not given, every base gets
# quality score 50)
sub fasta2fastq_hash{
    my $infile = shift;
    my $h_ref  = shift;

    my $id;
    my $seq;
    my $fake_q = 50;
    my $qual;          # to store fake qual scores

    if (scalar(@_) > 0){
        $fake_q = shift;
    }

    open F, $infile or croak "error opening $infile\n";
    my $s = join "", <F>;
    $s= substr($s,1);         # remove the first >
    $s =~ s/\n{2,}/\n/g;      # remove empty lines
    my @a = split /\n\>/, $s; # each element is id\nseq, where seq may have \n's
    undef $s;                 # free up some memory
    close F;

    # fill the hash
    foreach (@a){
        if (/(.*?\n)(.*)/s){
            $id = $1;
            chomp $id;

            $seq = $2;
            $seq =~ s/\n//g;

            # add the id/seq to hash
            $h_ref->{$id} = ();
            push @{$h_ref->{$id}}, ("$seq");

            # make fake quality scores and add them to hash
            $qual = "";

            for (my $i = 0; $i < length($seq); $i++){
                $qual .= mh12_utils::phred2fastq($fake_q);
            }

            push @{$h_ref->{$id}}, ("$qual");
        }
    }
}



#______________________________________________________________________________


# takes a fastq file and hash ref.  Fills the hash with id=>(array with
# 0 element: bases, 1st element: qual scores.
#
# NOTE: ASSUMES NO UNECESSARY LINE BREAKS IN THE FASTQ FILE: i.e. there
# are no \n's within the sequence or quality score strings (which is probably
# the spec for a fastq file, in any case)
#
# usage: fastq2length(fastq filename, hash reference, 1 | 0)
sub fastq2hash{
    my $infile = shift;
    my $h_ref  = shift;
    my $chop_header = shift;

    my $id;
    my $seq;

    open F, "$infile" or croak "error opening $infile\n";
    my $s = join "", <F>;
    $s= substr($s,1);         # remove the first @
    $s =~ s/\n{2,}/\n/g;      # remove empty lines
    my @a = split /\n\@/, $s; # each element is id\nseq\nid\nqual
    undef $s;                 # free up some memory
    close F;

    # fill the hash
    foreach (@a){
        if ( ($chop_header && /(.*?)\s+.*\n(.*)\n(.*)\n(.*)/)   ||  /(.*)\n(.*)\n(.*)\n(.*)/){
            $h_ref->{$1} = ();
            push @{$h_ref->{$1}}, ("$2", "$4");
        }
    }
}
#______________________________________________________________________________


# takes a fasta file and hash ref.  Fills the hash with id=>length
# usage: fasta2length(fasta filename, hash reference)
sub fasta2length{
    my $infile = shift;
    my $h_ref  = shift;

    my $id;
    my $seq;

    open F, $infile or croak "error opening $infile\n";
    my $s = join "", <F>;
    $s= substr($s,1);         # remove the first >
    $s =~ s/\n{2,}/\n/g;      # remove empty lines
    my @a = split /\n\>/, $s; # each element is id\nseq, where seq may have \n's
    undef $s;
    close F;

    # fill the hash
    foreach (@a){
        if (/(.*?\n)(.*)/s){
            $id = $1;
            chomp $id;

            $seq = $2;
            $seq =~ s/\n//g;

            $h_ref->{$id} = length $seq;
        }
    }
}
#______________________________________________________________________________


# usage: fasta2shortHeader(infile, outfile)
sub fasta2shortHeader{
    my $infile  = shift;
    my $outfile = shift;

    my %h;
    fasta2hash($infile, \%h);

    foreach(keys %h){
        if (/(.*?)\s.*/){
            $h{$1} = delete $h{$_};
        }
    }

    open F, ">$outfile" or die "Error opening $outfile\n";

    foreach (keys %h){
        print F ">$_\n$h{$_}\n";
    }

    close F;
}


# usage: fasta2singleLine(infile, outfile) or (infile, outfile, 1)
#
# takes infile, removes any unnecessary linebreaks and writes the result
# to outfile.  Using the 1 as a third argument will replace each newline
# with a space: useful for use on a qual file
sub fasta2singleLine {
    my $infile       = shift;
    my $outfile      = shift;
    my $insert_space = shift;

    # put whole file into an array, 1 sequence per element
    open F, $infile or die ("Error opening $infile");
    my $whole_file = join ("", <F>);
    close F;
    chomp $whole_file;
    $whole_file = substr ($whole_file,1);  # removes the > at the start of the file
    my @sequences = split("\n>", $whole_file);


    # write the output file
    open F, ">$outfile" or die ("Error opening $outfile");

    foreach (@sequences) {
        my @seq = split "\n";        # seq[0]=header, seq[1],seq[2],...=sequence
        print F ">$seq[0]\n";        # print seq header to file
        shift @seq;                  # take the header out of the array


        # put the sequence into 1 string (no newlines)
        my $out;

        if ($insert_space){
            $out = join " ", @seq;
        }
        else{
            $out = join "", @seq;
        }

        print F "$out\n";            # print sequence to file
    }

    close F;
}
#______________________________________________________________________________



# takes a fastq file and hash ref.  Fills the hash with id=>length
#
# NOTE: ASSUMES NO UNECESSARY LINE BREAKS IN THE FASTQ FILE: i.e. there
# are no \n's within the sequence or quality score strings (which is probably
# the spec for a fastq file, in any case)
#
# usage: fastq2length(fastq filename, hash reference)
sub fastq2length{
    my $infile = shift;
    my $h_ref  = shift;

    my $id;
    my $seq;

    open F, $infile or croak "error opening $infile\n";
    my $s = join "", <F>;
    $s= substr($s,1);         # remove the first @
    $s =~ s/\n{2,}/\n/g;      # remove empty lines
    my @a = split /\n\@/, $s; # each element is id\nseq\nid\nqual
    undef $s;                 # free up some memory
    close F;

    # fill the hash
    foreach (@a){
        if (/(.*)\n(.*)\n(.*)\n(.*)/){
            $h_ref->{$1} = length $2;
        }
    }
}
#______________________________________________________________________________


# puts the contents of a file into a hash - one line (chomped) per hash
# key, with key => 1, or key => nextline

# arg0: filename
# arg1: hash reference
# arg2: 0 for line=> 1;  1 for line => nextline
# arg3: 1 for remove 1st char from keys (useful for getting reads from fastas,
#       0 to not do so
sub file2hash {
    my $filename    = shift;
    my $hash_ref    = shift;
    my $next        = shift;
    my $rm_1st_char = shift;

    my $nextline;

    open F, $filename or croak "error opening file $filename in subroutine file2hash\n";
    while (my $line = <F>) {
        if ($line eq "\n"){next}  # skip empty lines
        chomp $line;
        if ($rm_1st_char){$line = substr $line, 1}

        if ($next){
            $nextline = <F>;
            while ($nextline eq "\n"){$nextline = <F>}  # skip empty lines
            chomp $nextline;
            $hash_ref->{$line} = $nextline;
        }
        else{
            $hash_ref->{$line} = 1;
        }
    }
}
#______________________________________________________________________________



# usage: phred2fastq(phred_score)
# returns the corresponding fastq score of phred_score
sub phred2fastq{
    my $p = shift;
    my $q = chr(($p<=93? $p : 93) + 33);
    return $q;
}
#______________________________________________________________________________


# usage: read2F_or_R(readname)
# returns: F or R
sub read2F_or_R{
    my $read = shift;

    if    ($read =~ /\.p/) {return "F"}
    elsif ($read =~ /\.q/) {return "R"}
    else{
        croak "error in subroutine read2F_or_R.  Readname $read not recognised\n";
    }
}
#______________________________________________________________________________


# usage: read2mate(readname)
# returns: name of the read's mate
sub read2mate {
    my $read = shift;

    if    ($read =~ /\.p/){$read =~ s/\.p/\.q/}
    elsif ($read =~ /\.q/){$read =~ s/\.q/\.p/}
    else{
        croak "error in subroutine read2mate. Readname $read not recognised\n";
    }
    return $read;
}
#______________________________________________________________________________


# usage: read2plate(readname)
# returns the plate name
sub read2plate {
    my $read  = shift;
    if ($read =~ /(.*)[a-z][0-9]+\..+/){
        return $1;
    }
    else{
        croak "error in subroutine read2plate.  Readname $read not recognised\n";
    }
}
#______________________________________________________________________________


# usage: read2template(readname)
# returns the template name
sub read2template {
    my $read  = shift;
    if ($read =~ /(.*)\..*/){
        return $1;
    }
    else{
        croak "error in subroutine read2template.  Readname $read not recognised\n";
    }
}
#______________________________________________________________________________


# usage: read2well(readname)
# returns the well ID
sub read2well {
    my $read  = shift;
    if ($read =~ /.*([a-z][0-9]+)\..+/){
        return $1;
    }
    else{
        croak "error in subroutine read2well.  Readname $read not recognised\n";
    }
}
#______________________________________________________________________________


# usage: reverse_complement(reference to string)
# returns the string reverse complemented
sub reverse_complement {
    my $s = shift;
    my $out = reverse $$s;

    $out =~ tr/acgtACGT/tgcaTGCA/;

    return $out;
}
#______________________________________________________________________________


# usage: subseq(string, start, end)
# returns substring from start to end inclusive.  If start > end, returns
# the reverse complement
# INDEX OF FIRST CHARACTER IS 1, NOT 0
sub subseq{
     my $s = shift;
     my $x = shift;
     my $y = shift;

    if ($x <= $y){
        $s = substr $s, $x-1, ($y-$x+1);
    }
    else{
        $s = substr $s, $y-1, ($x-$y+1);
        $s = mh12_utils::reverse_complement(\$s);
    }

    return $s;
}



#______________________________________________________________________________


# usage: mh12substr(string, start, end)
# returns substring from start to end inclusive.  If start > end, returns
# the characters in reverse order
# INDEX OF FIRST CHARACTER IS 1, NOT 0
sub mh12substr{
     my $s = shift;
     my $x = shift;
     my $y = shift;

    if ($x <= $y){
        $s = substr $s, $x-1, ($y-$x+1);
    }
    else{
        $s = reverse substr $s, $y-1, ($x-$y+1);
    }

    return $s;
}
#______________________________________________________________________________





1;

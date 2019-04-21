#!/usr/bin/perl

# selftest-compareresults.pl <new data> <reference data>

use strict;

my $retstatus = 0;

my $epsilon = 1.0E-10;
my $maxdelta = 0;

open my $local_fh, '<', $ARGV[0] or die "Can't open $ARGV[0]: $!";
my $newResults = resultsHash($local_fh);
close $local_fh;

open my $local_reffh, '<', $ARGV[1] or die "Can't open $ARGV[1]: $!";
my $refResults = resultsHash($local_reffh);
close $local_reffh;

foreach my $refKey (keys %$refResults) {
   if(exists($newResults->{$refKey})) {
      my $maxval = (abs($refResults->{$refKey})>=abs($newResults->{$refKey})?abs($refResults->{$refKey}):abs($newResults->{$refKey}));

      my $absdelta = abs($refResults->{$refKey} - $newResults->{$refKey});

      my $reldelta;
      my $delta;

      #print "item: $refKey\n";

      if($maxval > 0.0) {
      	$reldelta = $absdelta / $maxval;
      } else {
        $reldelta = 0.0;
      }

      $delta = ($absdelta <= $reldelta) ? $absdelta : $reldelta;
      if($delta > $maxdelta) {
         $maxdelta = $delta;
      }

      if($delta >= $epsilon) {
         print "$ARGV[0]: Significant difference for $refKey (reference: $refResults->{$refKey} new: $newResults->{$refKey} delta: $delta)\n";
         $retstatus = 1;
      } 
      delete $newResults->{$refKey};
   } else {
      print "$ARGV[0]: No corresponding value for $refKey\n";
      $retstatus = 1;
   }
}

foreach my $newKey (keys %$newResults) {
   print "$ARGV[0]: Found extra data item: $newKey -> $newResults->{$newKey}\n";
}

print "Maximum delta: $maxdelta\n";
exit $retstatus;

sub resultsHash {
   my $fd = shift;
   my %input;

   while( my $line = <$fd> ) {
      if ( $line =~ /(Estimate:|Assessed)/i ) {
         if($line =~ /^([^=]+) = ([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)$/ ) {
            my $label = $1;
            my $value = $2;

            $label =~ s/\s*$//;

            $input{"$label"} = $value;
            #print "$label -> $value\n";
         } 
      } 
   }

  return \%input; # return hash reference
}


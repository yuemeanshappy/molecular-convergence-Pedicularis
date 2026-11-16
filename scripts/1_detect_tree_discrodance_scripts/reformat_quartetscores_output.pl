#!/usr/bin/perl -w
use strict;

my $in = $ARGV[0];

my $tree;
open IN, "< $in" or die "Can't open $in!\n";
while (<IN>)	{
	s/\s//g;
	$tree .= $_;
}
close (IN);

while ($tree =~ /\[(.*?)\]/)	{
	my $ic = $1;
	my @ic = split /\;/, $ic;
	map { $ic[$_] =~ s/.*\:// }	(0..$#ic);
	$ic = join "/", @ic;
	$tree =~ s/\)([^\)]+)\[.*?\]/\)$ic\/$1/;
}

print "$tree\n";

exit;
#!/usr/bin/perl

# program to migrate power and temperature trace data for a 
# floorplan converted by hsconvert-flp.pl

# use warnings;

sub usage {
    print("usage: hsconvert-data.pl -o <original flp> -c <converted flp> -p <power file> (or)\n");
    print("       hsconvert-data.pl -o <original flp> -c <converted flp> -t <temp file>\n");
	print("migrates power/temperature trace data for a floorplan converted by hsconvert-flp.pl\n");
	print("prints output to stdout\n");
    print("-o <original flp>  -- floorplan file used as input to hsconvert-flp.pl (eg: ev6.flp.orig)\n");
    print("-c <converted flp> -- floorplan file containing the output of hsconvert-flp.pl (eg: ev6.flp)\n");
    print("-p <power file>    -- input power trace file corresponding to <original flp> (eg: gcc.ptrace.orig)\n");
    print("-t <temp file>     -- input temperature trace file corresponding to <original flp>\n");
    exit(1);
}

sub readnames {
	my $file = shift;
	my @strs;
	my @names;
	my $tail = 0;

	open(FILE, "<$file") || die("error: file $file could not be opened for reading\n");
	while (<FILE>) {
		# skip comments and empty lines
		next if (/^\s*#|^\s*$/);
		chomp;
		@strs = split(/\s+/);
		die ("error: wrong floorplan input format\n") if (@strs != 3 && @strs != 5);
		# skip connectivity information
		next if (@strs == 3);
		$names[$tail++] = $strs[0];
	}
	close(FILE);
	return @names;
}

# the no. of sub-blocks into which each block (that has been chopped) has been chopped
sub getchopcount {
	my $ofile = shift;
	my $cfile = shift;
	# original names
	my @onames = readnames($ofile);
	# converted names
	my @cnames = readnames($cfile);
	my %chopcount;
	my $i;
	my $j;

	for($i=0; $i < @onames; $i++) {
		for($j=0; $j < @cnames; $j++) {
			if($cnames[$j] =~ /^$onames[$i]_\d+$/) {
				if ($chopcount{$onames[$i]}) {
					$chopcount{$onames[$i]} += 1;
				} else {
					$chopcount{$onames[$i]} = 1;
				}
			}
		}
	}

	return %chopcount;
}

# main computation
usage() if (@ARGV != 6 || $ARGV[0] ne "-o" || ! -f $ARGV[1]  ||
			$ARGV[2] ne "-c" || ! -f $ARGV[3] || 
			$ARGV[4] ne "-p" && $ARGV[4] ne "-t" || ! -f $ARGV[5]);

%chopcount = getchopcount($ARGV[1], $ARGV[3]);

open (FILE, "<$ARGV[5]") || die("error: file $ARGV[5] could not be opened for reading\n");

# find the first non-empty line
while (<FILE>) {
	if (!/^\s*$/) {
		last;
	}
}
chomp;
@names = split(/\s+/);
die("error: empty trace file\n") if (!@names);

# the first line
for($i = 0; $i < @names-1; $i++) {
	if ($count = $chopcount{$names[$i]}) {
		for($j = 0; $j < $count; $j++) {
			print $names[$i]."_".$j."\t";
		}
	} else {
		print $names[$i]."\t";
	}
}
if ($count = $chopcount{$names[$i]}) {
	for($j = 0; $j < $count-1; $j++) {
		print $names[$i]."_".$j."\t";
	}
	print $names[$i]."_".$j."\n";
} else {
	print $names[$i]."\n";
}

# remaining lines
while(<FILE>) {
	# skip comments and empty lines
	next if (/^\s*#|^\s*$/);
	chomp;
	@vals = split(/\s+/);

	for($i = 0; $i < @vals-1; $i++) {
		if ($count = $chopcount{$names[$i]}) {
			# power values - divide proportionately
			if ($ARGV[4] eq "-p") {
				$vals[$i] /= $count;
			}
			for($j = 0; $j < $count; $j++) {
				print $vals[$i]."\t";
			}
		} else {
			print $vals[$i]."\t";
		}
	}
	if ($count = $chopcount{$names[$i]}) {
		if ($ARGV[4] eq "-p") {
			$vals[$i] /= $count;
		}
		for($j = 0; $j < $count-1; $j++) {
			print $vals[$i]."\t";
		}
		print $vals[$i]."\n";
	} else {
		print $vals[$i]."\n";
	}
}
close(FILE);


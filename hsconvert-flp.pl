#!/usr/bin/perl

# program to sanitize a floorplan with high aspect ratio blocks
# by chopping up of those blocks

# use warnings;

#params
$threshold = 1.0;

sub usage {
    print("usage: hsconvert-flp.pl [-a <aspect ratio>] <file>\n");
	print("sub-divides the high-aspect-ratio blocks and prints the resultant floorplan to stdout\n");
    print("[-a <aspect ratio>] -- approx aspect ratio threshold above which to sub-divide blocks (default 1.0)\n");
    print("<file>              -- input floorplan file (eg: ev6.flp.orig)\n");
    exit(1);
}

usage() if (@ARGV > 3 || !@ARGV % 2 || ! -f $ARGV[@ARGV-1]);

for($i=0; $i < @ARGV-1; $i+=2) {
	if ($ARGV[$i] eq "-a") {
		$threshold=$ARGV[$i+1];
		($threshold >= 1) || die("error: aspect ratio should be >= 1\n");
		next;
	}

	usage();
}

open (FILE, "<$ARGV[@ARGV-1]") || die("error: file $ARGV[@ARGV-1] could not be opened for reading\n");
while (<FILE>) {
	# skip comments and empty lines
	if (/^\s*#|^\s*$/) {
		print;
		next;
	}

	chomp;
	@strs = split(/\s+/);
	die ("error: wrong floorplan input format\n") if (@strs != 3 && @strs != 5);

	# block dimensions
	if (@strs == 5) {
		($name, $width, $height, $leftx, $bottomy) = @strs;
		$aspect = $height / $width;
		$vertical = 1;
		if ($aspect < 1) {
			$aspect = 1.0 / $aspect;
			$vertical = 0;
		}

		# block sub-division
		if ($aspect > $threshold) {
			$factor = int($aspect/$threshold + 0.5);
			# only one block after sub-division
			if ($factor <= 1) {
				print join("\t", @strs)."\n";
				next;
			}	
			$chopcount{$name} = $factor;
			if ($vertical) {
				$height /= $factor;
			} else {
				$width /= $factor;
			}
			for($i=0; $i < $factor; $i++) {
				printf($name."_".$i."\t%.6f\t%.6f\t%.6f\t%.6f\n", $width, $height, $leftx, $bottomy);
				if($vertical) {
					$bottomy += $height;
				} else {
					$leftx += $width;
				}	
			}
		# no need to sub-divide, print as-is
		} else {
			print join("\t", @strs)."\n";
		}
	# connectivity information
	} elsif(@strs == 3) {
		($name1, $name2, $density) = @strs;
		$count1 = $chopcount{$name1};
		$count2 = $chopcount{$name2};
		# update connectivity for sub-divided blocks
		if ($count1 && $count2) {
			for($i=0; $i < $count1; $i++) {
				for($j=0; $j < $count2; $j++) {
					printf($name1."_".$i."\t".$name2."_".$j."\t%.3f\n", $density);
				}
			}
		} elsif($count1) {
			for($i=0; $i < $count1; $i++) {
				printf($name1."_".$i."\t".$name2."\t%.3f\n", $density);
			}
		} elsif ($count2) {
			for($i=0; $i < $count2; $i++) {
				printf($name1."\t".$name2."_".$i."\t%.3f\n", $density);
			}	
		# print as-is
		} else {
			print join("\t", @strs)."\n";
		}
	}
}
close(FILE);


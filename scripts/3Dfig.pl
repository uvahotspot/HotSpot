#!/usr/bin/perl

# program to convert the floorplan file and FIG-like output 
# of HotSpot into a FIG file

#use warnings;

#params
#$maxx = 12600;
#$maxy = 9600;
$maxx = 9600;
$maxy = 12600;
$res = 1200;
#aspect ratio of blocks to give vertical labels
$skinny = 3;

sub usage () {
    print("usage: 3Dfig.pl [-a <area ratio>] [-f <fontsize>] [-s <nskip>] <file>\n");
    print("Reads in a layer configuration file and creates a .FIG file for each floorplan listed.\n");
    print("Requires that all referenced files are in the same directory as 3Dfig.pl.\n");
    print("[-a <area ratio>] -- approx page occupancy by figure (default 0.95)\n");
    print("[-f <fontsize>]   -- font size to be used (default 10)\n");
    print("[-s <nskip>]      -- no. of entries to be skipped in input (default 0)\n");
    print("<file>            -- input .lcf file (eg: example.lcf)\n");
    exit(1);
}

&usage() if (@ARGV > 7 || !@ARGV % 2 || ! -f $ARGV[@ARGV-1]);

$occupancy = 0.95;
$fontsize = 10;
$nskip = 0;

for($i=0; $i < @ARGV-1; $i+=2) {
    if ($ARGV[$i] eq "-a") {
        $occupancy=$ARGV[$i+1];
        next;
    }

    if ($ARGV[$i] eq "-f") {
        $fontsize=$ARGV[$i+1];
        next;
    }

    if ($ARGV[$i] eq "-s") {
        $nskip=$ARGV[$i+1];
        next;
    }

    &usage();
    exit(1);
}

open (LCFFILE, "<$ARGV[@ARGV-1]") || die("error: file $ARGV[@ARGV-1] could not be opened\n");
while(<LCFFILE>){
    next if (/^\s*#|^\s*$/);
    my($line) = $_;
    chomp($line);
    # Remove ^m DOS ending #
    $str =~ s/\r//g;
    next unless ($line =~ m/\.flp$/);
    push(@infile,$line);
}
close(LCFFILE);
$numfiles = @infile;


for($i = 0;$i<$numfiles;$i++){
    $fileout = substr($infile[$i],0,-3);
    $fileout = "$fileout"."fig";
    push(@outfiles,$fileout);
#	print "$fileout\n";
}


for($k=0;$k<$numfiles;$k++){
    print $outfiles[$k] . "\n";
#	print $infile[$k] . "\n";
#	print $numfiles . "\n";

    open (OUTFILE,">$outfiles[$k]")|| die("error: file $outfiles[$k] could not be opened\n");
    open (FILE, "<$infile[$k]") || die("error: file $infile[$k] could not be opened\n");
    #open (OUTFILE,">$fileout");
    #open (FILE, "<$ARGV[@ARGV-1]") || die("error: file $ARGV[@ARGV-1] could not be opened\n");
    $maxfig = -inf;
    $minfig = inf;
    $figinput = 0;
    while (<FILE>) {
        if (/FIG starts/) {
            $figinput = 1;
            last;
        }
    }

    # This is a HotSpot floorplan file and not an output of print_flp_fig. 
    # So let us generate the print_flp_fig output ourselves and save it in a 
    # temporary file.
    if (!$figinput) {
        seek(FILE, 0, 0);
        $timestamp = time();
        $file = "$infile[$i].$timestamp";
        open(NEWFILE, ">$file") || die("error: file $file could not be created\n");
        print(NEWFILE "FIG starts\n");
        while (<FILE>) {
            # skip comments and empty lines
            next if (/^\s*#|^\s*$/);
            chomp;
            @strs = split(/\s+/);
            $size = @strs;
            if (@strs != 3 && @strs != 5 && @strs !=7) {
                unlink($file);
                die ("error: wrong floorplan input format\n");
            }	
            # skip connectivity information
            next if (@strs == 3);
            ($name, $width, $height, $leftx, $bottomy, $sheat, $resist) = @strs;
            $rightx = $leftx + $width;
            $topy = $bottomy + $height;
            $size = @strs;
            if($size==7){
                printf(NEWFILE "%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", 
                    $leftx, $bottomy, $leftx, $topy, $rightx, $topy, $rightx, $bottomy, 
                    $leftx, $bottomy, $sheat, $resist);
            } else {
                printf(NEWFILE "%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", 
                    $leftx, $bottomy, $leftx, $topy, $rightx, $topy, $rightx, $bottomy, 
                    $leftx, $bottomy);			
            }

            printf(NEWFILE "%s\n", $name);
        }
        print(NEWFILE "FIG ends\n");
        close(NEWFILE);
        close(FILE);
        open(FILE, "<$file") || die("error: file $file could not be opened\n");
        while (<FILE>) {
            last if (/FIG starts/);
        }	
    }

    $pos=tell(FILE);
    $j=0; 
    while (<FILE>) {
        last if (/FIG ends/);
        next if (/[a-zA-Z_]|^\s*$/);
        if ($j < $nskip) {
            $j++;
            next;
        }
        chomp;
        @nums = split(/\s+/);
        $numsize = @nums;
        #	print "@nums Size:$numsize\n";
        if($numsize==12){
            $tempResist = pop(@nums);	
            $tempSheat = pop(@nums);	
        }

        for ($i=0; $i < @nums; $i++) {
            $maxfig = $nums[$i] if ($nums[$i] > $maxfig);
            $minfig = $nums[$i] if ($nums[$i] < $minfig);
        }
        if($numsize==12){
            push(@nums,$tempSheat);	
            push(@nums,$tempResist);	
        }

    }

    $maxfig-=$minfig;
    $scale = (($maxx < $maxy)?$maxx:$maxy) / $maxfig * sqrt($occupancy);
    $xorig = ($maxx-$maxfig*$scale)/2;
    $yorig = ($maxy-$maxfig*$scale)/2;

    print OUTFILE "#FIG 3.1\nPortrait\nCenter\nInches\n$res 2\n";
    seek(FILE, $pos, 0);

    $j=0;
    while (<FILE>) {
        last if (/FIG ends/);
        if ($j < $nskip * 2) {
            $j++;
            next;
        }	
        chomp;
        @coords = split(/\s+/);
        $size = @coords;
        if($size==12){
            $resist = pop(@coords);
            $sheat=pop(@coords);
        }
        next if ($#coords == -1);
        @coords = map($_-$minfig, @coords);
        @coords = map($_*$scale, @coords);
        $leftx = $rightx = $coords[0];
        $bottomy = $topy = $coords[1];
        for($i=2; $i < @coords; $i++) {
            if ($i % 2) {
                $bottomy = $coords[$i] if ($coords[$i] < $bottomy);
                $topy = $coords[$i] if ($coords[$i] > $topy);
            } else {
                $leftx = $coords[$i] if ($coords[$i] < $leftx);
                $rightx = $coords[$i] if ($coords[$i] > $rightx);
            }
        }
        for ($i=0; $i < @coords; $i++) {
            if ($i % 2) {
                $coords[$i] = int($maxy - $coords[$i] - $yorig);
            }	
            else {
                $coords[$i] = int($coords[$i] + $xorig);
            }	
        }
        printf OUTFILE ("2 2 0 1 -1 7 0 0 -1 0.000 0 0 0 0 0 %d\n", int (@coords)/2);
        print OUTFILE "\t@coords\n";
        $xpos = int($xorig + ($leftx + $rightx) / 2.0);
        $ypos = int($maxy - ($bottomy + $topy) / 2.0 - $yorig);
        $angle = 0;
        $angle = 1.5708 if (($topy - $bottomy) > $skinny * ($rightx - $leftx));
        $name = <FILE>;
        chomp($name);
        print OUTFILE "4 1 -1 0 0 0 $fontsize $angle 4 -1 -1 $xpos $ypos $name\\001\n";
        if($resist!=0 && $sheat!=0){
            $ypos = $ypos + $fontsize + 100;	
            print OUTFILE "4 1 -1 0 0 0 $fontsize $angle 4 -1 -1 $xpos $ypos ";
            printf OUTFILE ("Resistivity:%.3f\\001\n",$resist);
            $ypos = $ypos + 2*$fontsize + 100;	
            print OUTFILE "4 1 -1 0 0 0 $fontsize $angle 4 -1 -1 $xpos $ypos ";
            printf OUTFILE ("Specific Heat:%.3f\\001\n",$sheat);
        }
        $resist = 0;
        $sheat = 0;
    }
    close(FILE);
    close(OUTFILE);


    # delete the temporary file created
    unlink ($file) if (!$figinput);
}

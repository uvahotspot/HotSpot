#!/usr/bin/perl -w
#This script helps to generate an SVG figure.
#(More informaiton about SVG can be found at http://www.w3.org/TR/SVG/). 
#SVG figures can be zoomed without losing resolutions.
#
#USAGE: grid_thermal_map.pl <flp_file> <grid_temp_file> > <filename>.svg
#
#Use IE or SVG Viewer to open "<filename>.svg". May need to enable XML
#in your IE. Also, in Linux, ImageMagick 'convert' could be used to convert
#it to other file formats 
#(eg `convert -font Helvetica <filename>.svg <filename>.pdf`).
#
#Inputs: 
#	(1) A floorplan file with the same format as the other HotSpot .flp files.
#	(2) A list of grid temperatures, with the format of each line: 
#      <grid_number>	<grid_temperature>
#	Example floorplan file (example.flp) and the corresponding grid 
#	temperature file (example.t) are inlcuded in this release. The resulting
#	SVG figure (example.svg) is also included. 
#
#Acknowledgement: The HotSpot developers would like to thank Joshua Rosenbluh 
#at Grinnell College, Iowa, for his information and help on the SVG figures.

use strict;

sub usage () {
    print("usage: grid_thermal_map.pl <flp_file> <grid_temp_file> > <filename>.svg (or)\n");
    print("       grid_thermal_map.pl <flp_file> <grid_temp_file> <rows> <cols> > <filename>.svg\n");
    print("       grid_thermal_map.pl <flp_file> <grid_temp_file> <rows> <cols> <min> <max> > <filename>.svg\n");
	print("prints an 'SVG' format visual thermal map to stdout\n");
    print("<flp_file>       -- path to the file containing the floorplan (eg: ev6.flp)\n");
    print("<grid_temp_file> -- path to the grid temperatures file (eg: sample.t)\n");
    print("<rows>           -- no. of rows in the grid (default 64)\n");
    print("<cols>           -- no. of columns in the grid (default 64)\n");
    print("<min>            -- min. temperature of the scale (defaults to min. from <grid_temp_file>)\n");
    print("<max>            -- max. temperature of the scale (defaults to max. from <grid_temp_file>)\n");
    exit(1);
}

&usage() if (@ARGV != 2 && @ARGV != 4 && @ARGV != 6 || ! -f $ARGV[0] || ! -f $ARGV[1]);

#constants used throughout the program
	my $num_levels=21;			#number of colors used
	my $max_rotate=200; 		#maximum hue rotation
	my $floor_map_path=$ARGV[0];#path to the file containing floorplan
	my $temp_map_path=$ARGV[1]; #path to the grid temperature file 
	my $stroke_coeff=10**-7;	#used to tune the stroke-width
	my $stroke_opacity=0;		#used to control the opacity of the floor plan
	my $smallest_shown=10000;	#fraction of the entire chip necessary to see macro
	my $zoom=10**6;
	my $in_minx=0; my $in_miny=0; my $in_maxx=0; my $in_maxy=0; 
	my $txt_offset=100;
	my $x_bound;


#variables used throughout the program
my $fp  ;#the SVG to draw the floorplan
my $tm  ;#the SVG to draw the thermal map
my $defs;#definitions

my $min_x=10**10;my $max_x=0;
my $min_y=10**10;my $max_y=0;
my $tot_x;my $tot_y;

#Specify grid row and column here
my $row=64; my $col=64;
if (@ARGV >= 4) {
	$row = $ARGV[2];
	$col = $ARGV[3];
}

# define my palette with 21 RGB colors, from red to green to blue 
my @palette;

	$palette[0]='255,0,0';
	$palette[1]='255,51,0';
	$palette[2]='255,102,0';
	$palette[3]='255,153,0';
	$palette[4]='255,204,0';
	$palette[5]='255,255,0';
	$palette[6]='204,255,0';
	$palette[7]='153,255,0';
	$palette[8]='102,255,0';
	$palette[9]='51,255,0';
	$palette[10]='0,255,0';
	$palette[11]='0,255,51';
	$palette[12]='0,255,102';
	$palette[13]='0,255,153';
	$palette[14]='0,255,204';
	$palette[15]='0,255,255';
	$palette[16]='0,204,255';
	$palette[17]='0,153,255';
	$palette[18]='0,102,255';
	$palette[19]='0,51,255';
	$palette[20]='0,0,255';
	
{#generate the SVG for the floorplan

	#declare variables to be used locally
	my @AoA_flp;#Array of arrays
	my @AoA_grid;
	my $min_t=1000;
	my $max_t=0;
	
	my $unit_cnt=0;
	
	#this section reads the floor map file
	{my $num='\s+([\d.]+)';my $dumb_num='\s+[\d.]+';
	open(FP, "< $floor_map_path")    or die "Couldn't open $floor_map_path for reading: $!\n";
	#input file order: instance name, width, height, minx, miny

	while (<FP>) {
			if (/^(\S+)$num$num$num$num/)
			{	
				$in_minx=$4*$zoom; $in_miny=$5*$zoom; $in_maxx=($4+$2)*$zoom; $in_maxy=($5+$3)*$zoom;
				$min_x=$in_minx if $in_minx<$min_x;$max_x=$in_maxx if $in_maxx>$max_x;
				$min_y=$in_miny if $in_miny<$min_y;$max_y=$in_maxy if $in_maxy>$max_y;
#				$min_t=$6 if $6<$min_t;$max_t=$6 if $6>$max_t;
				push @AoA_flp, [ ($1,$in_minx,$in_miny,$in_maxx,$in_maxy,$6) ] ;
			}								
	}
	
	close(FP);
	}
	$tot_x=$max_x-$min_x;
	$tot_y=$max_y-$min_y;
	$x_bound=$max_x*1.2;
	
	my $grid_h=$tot_y/$row;
	my $grid_w=$tot_x/$col;
	my ($grid_minx,$grid_miny,$grid_maxx,$grid_maxy);
	
	#this section reads the temperature map file
	{my $num1='\s+([\d.]+)';my $dumb_num1='\s+[\d.]+';
	open(TM, "< $temp_map_path")    or die "Couldn't open $temp_map_path for reading: $!\n";
	#input file order: grid#, temperature

	while (<TM>) {
		if (/^(\S+)$num1/)
		{
			$grid_minx=($1%$col)*$grid_w;
			$grid_maxx=($1%$col+1)*$grid_w;
			$grid_miny=(int($1/$col))*$grid_h;
			$grid_maxy=(int($1/$col)+1)*$grid_h;
			$min_t=$2 if $2<$min_t;$max_t=$2 if $2>$max_t;
			push @AoA_grid, [ ($1,$grid_minx,$grid_miny,$grid_maxx,$grid_maxy,$2) ] ;
		}								
	}
	close(TM);
	}
	
# if upper and lower limits need to be specified, do it here
if (@ARGV == 6) {
	$min_t=$ARGV[4];
	$max_t=$ARGV[5];
}
	
	#draw the grid temperatures
	$tm='<g id="floorplan" style="stroke: none; fill: red;">'."\n";
	{
	my ($w1,$h1, $level);
	foreach (@AoA_grid){
		$w1=(@{$_}[3])-(@{$_}[1]);
		$h1=(@{$_}[4])-(@{$_}[2]);
		if ($w1>$tot_x/$smallest_shown && $h1>$tot_y/$smallest_shown){
			if ($max_t > $min_t) {
				$level=int(($max_t-(@{$_}[5]))/($max_t-$min_t)*($num_levels-1));
			} else {
				$level = 0;
			}
			$tm.="\t".'<rect x="'.@{$_}[1] .'" y="'. @{$_}[2] .
			'" width="'.$w1 .'" height="'.$h1 .
			'" style="fill:rgb(' .$palette[$level].')" />'."\n";
		}
	}
	}

	#draw the floorplan
	{
	my ($w,$h, $start_y, $end_y, $txt_start_x, $txt_start_y);
	foreach (@AoA_flp){
		$unit_cnt += 1;
		$w=(@{$_}[3])-(@{$_}[1]);
		$h=(@{$_}[4])-(@{$_}[2]);
		$start_y=$tot_y-@{$_}[2]-$h;
		$end_y=$tot_y-@{$_}[2];
		$txt_start_x=@{$_}[1]+$txt_offset;
		$txt_start_y=$start_y+2*$txt_offset;
		if ($w>$tot_x/$smallest_shown && $h>$tot_y/$smallest_shown){
			$fp.="\t".'<line x1="'.@{$_}[1] .'" y1="'. $start_y .
			'" x2="'. @{$_}[3] .'" y2="'. $start_y .
			'" style="stroke:black;stroke-width:30" />'."\n";
			
			$fp.="\t".'<line x1="'.@{$_}[1] .'" y1="'. $start_y .
			'" x2="'. @{$_}[1] .'" y2="'. $end_y .
			'" style="stroke:black;stroke-width:30" />'."\n";
			
			$fp.="\t".'<line x1="'.@{$_}[3] .'" y1="'. $start_y .
			'" x2="'. @{$_}[3] .'" y2="'. $end_y .
			'" style="stroke:black;stroke-width:30" />'."\n";
			
			$fp.="\t".'<line x1="'.@{$_}[1] .'" y1="'. $end_y .
			'" x2="'. @{$_}[3] .'" y2="'. $end_y .
			'" style="stroke:black;stroke-width:30" />'."\n";
			
			$fp.="\t".'<text x="'.$txt_start_x .'" y="'. $txt_start_y .
			'" fill="black" text_anchor="start" style="font-size:180" > '. @{$_}[0] .' </text>'."\n";
		}
	}
	}

# draw the color scale bar	
	{
	my $i;
	my $txt_ymin;
	my $w2=$max_x*0.05;
	my $h2=$max_y*0.025;
	my $clr_xmin=$max_x*1.1;
	my $clr_ymin=$max_y*0.05;
	my $scale_xmin=$max_x*1.05;
	my $scale_value;
	my $final_scale_value=$min_t+(1/$num_levels)*($max_t-$min_t);
	$final_scale_value=~s/^(\d+)\.(\d)(\d)(\d)\d+/$1\.$2$3$4/;
	
	for ($i=0; $i<$num_levels; $i++) {
		if ($w2>$tot_x/$smallest_shown && $h2>$tot_y/$smallest_shown){
#			$level=int(($max_t-(@{$_}[5]))/($max_t-$min_t)*($num_levels-1));
			$fp.="\t".'<rect x="'.$clr_xmin .'" y="'. $clr_ymin .
			'" width="'.$w2 .'" height="'.$h2 .
			'" style="fill:rgb(' .$palette[$i].'); stroke:none" />'."\n";
			if ($i%3==0) {
				$txt_ymin=$clr_ymin+$h2*0.5;
				$scale_value=($max_t-$min_t)*(1-$i/($num_levels-1))+$min_t;
				$scale_value=~s/^(\d+)\.(\d)(\d)\d+/$1\.$2$3/;
				$fp.="\t".'<text x="'.$scale_xmin .'" y="'. $txt_ymin.
				'" fill="black" text_anchor="start" style="font-size:250" > '. $scale_value .' </text>'."\n";
			}
		}
		$clr_ymin+=$h2;
	}
	$min_t=~s/^(\d+)\.(\d)(\d)\d+/$1\.$2$3/;
	$fp.="\t".'<text x="'.$scale_xmin .'" y="'. $clr_ymin .
	'" fill="black" text_anchor="start" style="font-size:250" > '. $min_t .' </text>'."\n";
	}
	
	$fp.="</g>\n";
}
#svg header and footer
my $svgheader= <<"SVG_HEADER";
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN"
    "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg width="900px" height="500px"
    viewBox="$min_x $min_y $x_bound $max_y">
<title>Sample Temperature Map For HotSpot Grid Model</title>
SVG_HEADER

my $svgfooter= <<"SVG_FOOTER";
</svg>
SVG_FOOTER

print $svgheader.$tm.$fp.$svgfooter;
#writes out the header, definitions, thermal map, floor plan and footer

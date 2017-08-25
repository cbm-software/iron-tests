#!/usr/bin/perl

#use strict;
use warnings;


my ($open_exnode_file, $filename, $force, $force2, $force3, );

# Auszuwertende Knoten
my @biceps_bottom_nodes=(1..33); # Right now meaningless, reads in every node

print "Evaluating nodes: @biceps_bottom_nodes \n";
print "#######################################\n";

# Zeitschritte
my @disps = (5..5);

open WRITE_RESULTS, ">displacements.dat" or die "$!"; 
printf WRITE_RESULTS "# time   force\n";

  foreach my $disp (@disps)
  {
    $filename = sprintf("Example_%d", $disp);
    $open_exnode_file = "./$filename.part0.exnode";
    print "Reading: $open_exnode_file \n" ;
    $force = call_biceps_bottom(@biceps_bottom_nodes);

    printf WRITE_RESULTS "%3.7f  %2.6e \n", ($disp), $force;
    
    $force2 = sprintf("my_displacement_%d", $disp);
    $force3 = "cp all_fields $force2.txt";
    system($force3);
    
  }

close WRITE;







#system("mv all_fields field_*");
#system("rm all_fields field_* reac_forces");

sub call_biceps_bottom{

    my @biceps_bottom_nodes=@_;
	
    my $write_file="all_fields";
    my $line_no=-1;
    my $node_no;
    my $first=0;
#    my $write4="field_4";
    my $sum_forces_x=0;
    my $CoV=0;
    my $ave;
    my $std;
    my $k=0;
    my $line;

# im prinzip nur die daten zu knoten aus den exnode files raussuchen und in die reaktionskräfte incl. ableitugen in all_fields sammeln 
    open READ, "$open_exnode_file" or die "$!"; #wurde von start_cmiss übergeben
    open WRITE, ">$write_file" or die "$!"; #all_fields
    while(my $line=<READ>){
      #m=The match operator; suche nach Node, \s - string \d - number      
	    if($line=~m/^ Node:[\s]+([\d]+)/) #wenn das zutrifft dann tue:
	      {$line_no=0;$node_no=$1;$first=1}
	    if(($line_no==4) && ($first==1)) # Read x-displacement (field 3)
	    {
	      $line = trim($line);
	      print WRITE $line,"\n";
	    }
	    $line_no++;
    }
    close READ;
    close WRITE;


    return $CoV;

} #ende subroutine


sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };


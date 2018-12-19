#!/usr/bin/perl -w

use strict; use warnings;
use Algorithm::Combinatorics qw(combinations);


if (@ARGV != 2) {
  die "\nUsage: GetNucFrequency.pl TableFile K-mer > Outfile\n\n";
}

my $file1=shift;
my $kmero=shift;
my @all_kmer=split /,/, $kmero;
my %khash=();


open(IN, "<$file1") or die ("Couldn't open file $file1\n");
my $first =<IN>;
chomp ($first);
my @names=split /\s+/, $first;
while (my $l=<IN>) {
    chomp($l);
    my @line = split /\s+/, $l;
    my $key=$line[0];
    #my @lettersplit=split // $key;
    for (my $s=1; $s<@line; $s++){
        $khash{$names[$s]}{$key}=$line[$s];
    }
}
close IN;

my @samples=sort keys %khash;


print "Matrix\t";
print join "\t", @samples;
print "\n";

foreach my $h (@all_kmer){
  my @kmer = split //, $h;
  my @poskmer=();
  for (my $i=0; $i<@kmer; $i++){
    $poskmer[$i]=$i;
  }
  print "$h";

  foreach my $l (@samples){
      #print "En muestra $l\n";
      my $ksize=scalar(@kmer);
      my $current=0;
      my $cumulative=1;
      unless ($khash{$l}{$h} && $khash{$l}{$h}>0) {
        $cumulative=0;
        print "\t$cumulative";
        next;
      }

      for (my $i=$ksize; $i>=1; $i--){
          my $iter = combinations(\@poskmer, $i);
          while (my $c = $iter->next) {
              my $testk="";
              my $startEnd = ($c->[-1]-$c->[0])+1;
              my $size = scalar(@{$c});
              #print "Startend es $startEnd y size $size\n";
              if ($startEnd > $size) {
                  #Para cada posible K-mero debe sumar las probabilidades y esa suma es el kmero por el cual va a multiplicar.
                  #For each possible K-mer you must add the probabilities and that sum is the cell by which you will multiply.
                  my $count=$c->[0];
                  my $t=0;
                  my $numN=0;
                  while ($t<@{$c}){
                    if ($c->[$t]==$count) {
                      $testk.=$kmer[$c->[$t]];
                      $count++;
                      $t++;
                    }else{
                      $testk.="N";
                      $count++;
                      $numN++;
                    }
                  }
                  # Ahora tiene que definir el valor del Freq con N basado en la suma de todos los componentes de N.
                  # Now you have to define the Freq value with N based on the sum of all the N components
                  unless ($khash{$l}{$testk}){
                    $khash{$l}{$testk}=0;
                    repN($testk, $testk, $l);
                  }
              }else{
                  for (my $t=$c->[0]; $t<=$c->[-1]; $t++){
                    $testk.=$kmer[$t];
                  }
              }
              if ($current==0) {
                  #print "Va a multiplicar por freq >$testk<\n";
                  $cumulative*=$khash{$l}{$testk};
              }else {
                  #print "Va a dividir por freq >$testk<\n";
                  #if ($khash{$l}{$testk}==0) {
                  #  $cumulative=0;
                  #}else{
                    $cumulative= $cumulative/=$khash{$l}{$testk};
                  #}
              }
          }
          $current= $current==0 ? 1 : 0;
      }
      print "\t$cumulative";
  }
  print "\n";
}


sub repN {
  my $n=shift(@_);
  my $original=shift(@_);
  my $sample=shift(@_);
  if ($n=~/N/) {
    my $a=$n;
    $a=~s/N/A/;
    if ($a=~/N/) {
      repN($a, $original, $sample);
    }else{
      $khash{$sample}{$original}+=$khash{$sample}{$a};
      #print "$a\t$original\n";

    }
    my $t=$n;
    $t=~s/N/T/;
    if ($t=~/N/) {
      repN($t, $original, $sample);
    }else{
      $khash{$sample}{$original}+=$khash{$sample}{$t};
      #print "$t\t$original\n";
    }
    my $c=$n;
    $c=~s/N/C/;
    if ($c=~/N/) {
      repN($c, $original, $sample);
    }else{
      $khash{$sample}{$original}+=$khash{$sample}{$c};
      #print "$c\t$original\n";
    }
    my $g=$n;
    $g=~s/N/G/;
    if ($g=~/N/) {
      repN($g, $original, $sample);
    }else{
      $khash{$sample}{$original}+=$khash{$sample}{$g};
      #print "$g\t$original\n";
    }
  }
}

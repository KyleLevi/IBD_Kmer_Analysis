#!/usr/bin/perl -w

use strict;

if (@ARGV != 1) {
  die "\nUsage: GetNucFrequency.pl FastaFile > Outfile\n\n";
}


my $file1=shift;
my @array=('AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC','AAA','AAT','AAG','AAC','ATA','ATT','ATG','ATC','AGA','AGT','AGG','AGC','ACA','ACT','ACG','ACC','TAA','TAT','TAG','TAC','TTA','TTT','TTG','TTC','TGA','TGT','TGG','TGC','TCA','TCT','TCG','TCC','GAA','GAT','GAG','GAC','GTA','GTT','GTG','GTC','GGA','GGT','GGG','GGC','GCA','GCT','GCG','GCC','CAA','CAT','CAG','CAC','CTA','CTT','CTG','CTC','CGA','CGT','CGG','CGC','CCA','CCT','CCG','CCC','AAAA','AAAT','AAAG','AAAC','AATA','AATT','AATG','AATC','AAGA','AAGT','AAGG','AAGC','AACA','AACT','AACG','AACC','ATAA','ATAT','ATAG','ATAC','ATTA','ATTT','ATTG','ATTC','ATGA','ATGT','ATGG','ATGC','ATCA','ATCT','ATCG','ATCC','AGAA','AGAT','AGAG','AGAC','AGTA','AGTT','AGTG','AGTC','AGGA','AGGT','AGGG','AGGC','AGCA','AGCT','AGCG','AGCC','ACAA','ACAT','ACAG','ACAC','ACTA','ACTT','ACTG','ACTC','ACGA','ACGT','ACGG','ACGC','ACCA','ACCT','ACCG','ACCC','TAAA','TAAT','TAAG','TAAC','TATA','TATT','TATG','TATC','TAGA','TAGT','TAGG','TAGC','TACA','TACT','TACG','TACC','TTAA','TTAT','TTAG','TTAC','TTTA','TTTT','TTTG','TTTC','TTGA','TTGT','TTGG','TTGC','TTCA','TTCT','TTCG','TTCC','TGAA','TGAT','TGAG','TGAC','TGTA','TGTT','TGTG','TGTC','TGGA','TGGT','TGGG','TGGC','TGCA','TGCT','TGCG','TGCC','TCAA','TCAT','TCAG','TCAC','TCTA','TCTT','TCTG','TCTC','TCGA','TCGT','TCGG','TCGC','TCCA','TCCT','TCCG','TCCC','GAAA','GAAT','GAAG','GAAC','GATA','GATT','GATG','GATC','GAGA','GAGT','GAGG','GAGC','GACA','GACT','GACG','GACC','GTAA','GTAT','GTAG','GTAC','GTTA','GTTT','GTTG','GTTC','GTGA','GTGT','GTGG','GTGC','GTCA','GTCT','GTCG','GTCC','GGAA','GGAT','GGAG','GGAC','GGTA','GGTT','GGTG','GGTC','GGGA','GGGT','GGGG','GGGC','GGCA','GGCT','GGCG','GGCC','GCAA','GCAT','GCAG','GCAC','GCTA','GCTT','GCTG','GCTC','GCGA','GCGT','GCGG','GCGC','GCCA','GCCT','GCCG','GCCC','CAAA','CAAT','CAAG','CAAC','CATA','CATT','CATG','CATC','CAGA','CAGT','CAGG','CAGC','CACA','CACT','CACG','CACC','CTAA','CTAT','CTAG','CTAC','CTTA','CTTT','CTTG','CTTC','CTGA','CTGT','CTGG','CTGC','CTCA','CTCT','CTCG','CTCC','CGAA','CGAT','CGAG','CGAC','CGTA','CGTT','CGTG','CGTC','CGGA','CGGT','CGGG','CGGC','CGCA','CGCT','CGCG','CGCC','CCAA','CCAT','CCAG','CCAC','CCTA','CCTT','CCTG','CCTC','CCGA','CCGT','CCGG','CCGC','CCCA','CCCT','CCCG','CCCC');

my $seq1="";
my %counts=();
my %counts_1=();
my %counts_2=();
my %counts_3=();
my %counts_4=();
my $total_mono=0;
my $total_bi=0;
my $total_tri=0;
my $total_tetra=0;

open (IN, "<$file1") or die ("Couldn't open file $file1\n");
my $name1=$file1;
chomp ($name1);

while (my $i=<IN>){
  next unless ($i =~ /\w+/);
  chomp($i);
  if ($i =~ /^>/){
    unless ($seq1 eq ""){
      $seq1 =~ s/[^ATCG]//g;
      &process_nuc($seq1, $name1);
    }
    $seq1="";
  }else{
    $seq1.=uc($i);
  }
}
close IN;
$seq1 =~ s/[^ATCG]//g;
&process_nuc($seq1, $name1);


print "Matrix_";
print scalar(@array);
for (my $k=0; $k<@array; $k++){
  print "\t$array[$k]";
}print "\n";

my %norm_1=();
my %norm_2=();
my %norm_3=();
my %norm_4=();

foreach my $k (keys (%counts)){
  print "$k";

  my $value=0;
  $norm_1{'A'}=$counts_1{$k}{'A'}/$total_mono;
  $norm_1{'T'}=$counts_1{$k}{'T'}/$total_mono;
  $norm_1{'C'}=$counts_1{$k}{'C'}/$total_mono;
  $norm_1{'G'}=$counts_1{$k}{'G'}/$total_mono;

  for (my $j=0; $j<@array; $j++){
    my $eval=$array[$j];
    my @ev=split //, $eval;  
    if (length($eval)==2){
      $norm_2{$ev[0]}{$ev[1]}=$counts_2{$k}{$ev[0]}{$ev[1]}/$total_bi if ($counts_2{$k}{$ev[0]}{$ev[1]});
      $norm_2{$ev[0]}{'X'}{$ev[1]}=$counts_2{$k}{$ev[0]}{'X'}{$ev[1]}/$total_tri if ($counts_2{$k}{$ev[0]}{'X'}{$ev[1]});
      $norm_2{$ev[0]}{'X'}{'X'}{$ev[1]}=$counts_2{$k}{$ev[0]}{'X'}{'X'}{$ev[1]}/$total_tetra if ($counts_2{$k}{$ev[0]}{'X'}{'X'}{$ev[1]});
    }
    if (length($eval)==3){
      $norm_3{$ev[0]}{$ev[1]}{$ev[2]}=$counts_3{$k}{$ev[0]}{$ev[1]}{$ev[2]}/$total_tri if ($counts_3{$k}{$ev[0]}{$ev[1]}{$ev[2]});
      $norm_3{$ev[0]}{'X'}{$ev[1]}{$ev[2]}=$counts_3{$k}{$ev[0]}{'X'}{$ev[1]}{$ev[2]}/$total_tetra if ($counts_3{$k}{$ev[0]}{'X'}{$ev[1]}{$ev[2]});
      $norm_3{$ev[0]}{$ev[1]}{'X'}{$ev[2]}=$counts_3{$k}{$ev[0]}{$ev[1]}{'X'}{$ev[2]}/$total_tetra if ($counts_3{$k}{$ev[0]}{$ev[1]}{'X'}{$ev[2]});
    }
    if (length($eval)==4){
      $norm_4{$ev[0]}{$ev[1]}{$ev[2]}{$ev[3]}=$counts_4{$k}{$ev[0]}{$ev[1]}{$ev[2]}{$ev[3]}/$total_tetra if ($counts_4{$k}{$ev[0]}{$ev[1]}{$ev[2]}{$ev[3]});
    }
  }

  #print "The total mono counts are $total_mono the total bi counts are $total_bi. Freq of A: $counts_1{$k}{'A'} => $norm_1{'A'} and T: $counts_1{$k}{'T'} => $norm_1{'T'}; so di nuc AA is $counts_2{$k}{'A'}{'A'} and norm is $norm_2{'A'}{'A'}; while AT is $counts_2{$k}{'A'}{'T'} and norm is $norm_2{'A'}{'T'}\n\n\n";

  for (my $j=0; $j<@array; $j++){
    my $eval=$array[$j];
    my @ev=split //, $eval;
    if (length($eval)==2){
      if ($norm_2{$ev[0]}{$ev[1]}){
	$value = sprintf("%.5f", ($norm_2{$ev[0]}{$ev[1]}/($norm_1{$ev[0]}*$norm_1{$ev[1]})));
      }else{ 
	$value = sprintf("%.5f", 0);
      }
    }elsif (length($eval)==3){
      if ($norm_3{$ev[0]}{$ev[1]}{$ev[2]}){
	$value = sprintf("%.5f", (($norm_3{$ev[0]}{$ev[1]}{$ev[2]}*$norm_1{$ev[0]}*$norm_1{$ev[1]}*$norm_1{$ev[2]})/($norm_2{$ev[0]}{$ev[1]}*$norm_2{$ev[1]}{$ev[2]}*$norm_2{$ev[0]}{'X'}{$ev[2]})));
      }else{ 
	$value = sprintf("%.5f", 0);
      }
    }elsif (length($eval)==4){
      if ($norm_4{$ev[0]}{$ev[1]}{$ev[2]}{$ev[3]}){
	$value = sprintf("%.5f", (($norm_4{$ev[0]}{$ev[1]}{$ev[2]}{$ev[3]}*$norm_2{$ev[0]}{$ev[1]}*$norm_2{$ev[0]}{'X'}{$ev[2]}*$norm_2{$ev[0]}{'X'}{'X'}{$ev[3]}*$norm_2{$ev[1]}{$ev[2]}*$norm_2{$ev[1]}{'X'}{$ev[3]}*$norm_2{$ev[2]}{$ev[3]})/($norm_3{$ev[0]}{$ev[1]}{$ev[2]}*$norm_3{$ev[1]}{$ev[2]}{$ev[3]}*$norm_3{$ev[0]}{'X'}{$ev[2]}{$ev[3]}*$norm_3{$ev[0]}{$ev[1]}{'X'}{$ev[3]}*$norm_1{$ev[0]}*$norm_1{$ev[1]}*$norm_1{$ev[2]}*$norm_1{$ev[3]})));
      }else{ 
	$value = sprintf("%.5f", 0);
      }
    }
    print "\t$value";
  }
  print "\n";
}


sub process_nuc {
  my $seq=$_[0];
  my $name= $_[1];
  my $revcomp="";
  my $mono=();
  my $bi="";
  my $tri="";
  my $tetra="";


  $total_mono+=length($seq)*2;
  $total_bi+=(length($seq)-1)*2;
  $total_tri+=(length($seq)-2)*2;
  $total_tetra+=(length($seq)-3)*2;
  $revcomp=reverse($seq);
  $revcomp=~ tr/ATCG/TAGC/;
  $counts{$name}=1;
  while ($seq){
    if (length($seq) > 3){
      $tetra=substr($seq, 0, 4);
      my @t = split //, $tetra;
      $counts_4{$name}{$t[0]}{$t[1]}{$t[2]}{$t[3]}++;
      $counts_2{$name}{$t[0]}{'X'}{'X'}{$t[3]}++;
      $counts_3{$name}{$t[0]}{'X'}{$t[2]}{$t[3]}++;
      $counts_3{$name}{$t[0]}{$t[1]}{'X'}{$t[3]}++;
      $tetra=substr($revcomp, 0, 4);
      @t = split //, $tetra;
      $counts_4{$name}{$t[0]}{$t[1]}{$t[2]}{$t[3]}++;
      $counts_2{$name}{$t[0]}{'X'}{'X'}{$t[3]}++;
      $counts_3{$name}{$t[0]}{'X'}{$t[2]}{$t[3]}++;
      $counts_3{$name}{$t[0]}{$t[1]}{'X'}{$t[3]}++;
    }
    if (length($seq) > 2){
      $tri=substr($seq, 0, 3);
      my @tr = split //, $tri;
      $counts_3{$name}{$tr[0]}{$tr[1]}{$tr[2]}++;
      $counts_2{$name}{$tr[0]}{'X'}{$tr[2]}++;
      $tri=substr($revcomp, 0, 3);
      @tr = split //, $tri;
      $counts_3{$name}{$tr[0]}{$tr[1]}{$tr[2]}++;
      $counts_2{$name}{$tr[0]}{'X'}{$tr[2]}++;
    }
    if (length($seq) > 1){
      $bi=substr($seq, 0, 2);
      my @b=split //, $bi; 
      $counts_2{$name}{$b[0]}{$b[1]}++;
      $bi=substr($revcomp, 0, 2);
      @b=split //, $bi; 
      $counts_2{$name}{$b[0]}{$b[1]}++;
    }
    $mono= substr $seq, 0, 1 , ""; 
    $counts_1{$name}{$mono}++;
    $mono= substr $revcomp, 0, 1 , ""; 
    $counts_1{$name}{$mono}++;
  }
  return 1;
}

#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Spec::Functions qw/rel2abs/;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Data::Dumper;

my $programe_dir=basename($0);
my $path=dirname($0);
my $ver    = "2.1";
my $Write_Date   = "2021-06-11";
my $BEGIN  = time();
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($snplist,$geno,$vcf, $outdir, $prefix ,$group, $if_nmfst ,$run, $size );
GetOptions(
			"h|?" =>\&help,
			"i:s" =>\$snplist,
			"t:s" =>\$geno,
			"v:s" =>\$vcf,
			"g:s" => \$group,
			"k:s" =>\$prefix,
			"r"   =>\$run,
			"od:s"=>\$outdir,
			"nm"  =>\$if_nmfst,
			"s:i" =>\$size
			) || &help;
&help unless ($snplist or $geno or $vcf);
&help unless ($outdir && $prefix);

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";
##########################################################
my $od = rel2abs($outdir);
system("mkdir -p $od") unless(-d "$od");

my %code1=(
"A"=>"AA","T"=>"TT","C"=>"CC","G"=>"GG",
"R"=>"AG","Y"=>"CT","M"=>"AC","K"=>"GT","S"=>"CG","W"=>"AT",
"H"=>"ACT","B"=>"CGT","V"=>"ACG","D"=>"AGT");
my %base = ('A C'=>'M','C A'=>'M','A G'=>'R','G A'=>'R','A T'=>'W','T A'=>'W',
'C G'=>'S','G C'=>'S','C T'=>'Y','T C'=>'Y','G T'=>'K','T G'=>'K',
'A A'=>'A','C C'=>'C','T T'=>'T','G G'=>'G');
my %Homcode = ("A"=>"AA","T"=>"TT","C"=>"CC","G"=>"GG",);
my %Hetcode = ("R"=>"AG","Y"=>"CT","M"=>"AC","K"=>"GT","S"=>"CG","W"=>"AT","?"=>"NN",
"H"=>"ACT","B"=>"CGT","V"=>"ACG","D"=>"AGT","N"=>"NN","-"=>"NN",);
my @bit = ("A", "T", "C", "G", "R", "Y", "M", "K", "S", "W");
my %na  = ("AA",1, "TT",1, "CC",1, "GG",1);

#print Dumper(\%na); die;
##########################################################

info("1.reading group file...");
my (%genotype,%loci, %num);
$group = (defined $group) ? $group : 1;
my %group = &deal_group($group);

my $flag  = defined $snplist ? "snp" : 
            defined $geno    ? "geno": "vcf";
my $infile= defined $snplist ? $snplist :
            defined $geno    ? $geno : $vcf; 

##########################################################
info("2.$flag trans to genotype file...");
&deal_geno($infile,$od,$prefix,\%group,$flag,\%genotype); 

my @index = qw/ae Ho He nei shi pic/;
if ($run){
	for my $g (keys %genotype){
		open my $out,">$outdir/loci_of_$prefix.$g.index.xls" or die $!;
		print $out join("\t","chr","loci","Number","maf",@index,"Na")."\n";
		$loci{$g} = $out;
	}
}

print Dumper(\%genotype); 

##########################################################
info("3.genetic index calculate...");

my %stat;
my %na_stat;
for my $g ( sort {$a cmp $b} keys %genotype){
	open IN, $genotype{$g} or die $!; 
	print "group: $g\n";
	while(<IN>){
		chomp;
		if ($.==1){
			next;
		}else{
			my ($Chr, $pos, @genotype) = split(/\t/, $_);  
			$stat{$g}{sample_count} = scalar @genotype;
			my $het = 0 ; my @notN;
			for my $type (@genotype){
				unless ($type eq "N"){
					$het++ if (exists $Hetcode{$type});
					push(@notN, $code1{$type});
				}
			}
			my $loc = "$Chr\t$pos";
			&deal_allele(\@notN,$g,$het,$loc) if (scalar @notN >1);
			map{$na_stat{$g}{na} += $na{$_} ? $na{$_} : 2; $na_stat{$g}{num} += 1}@notN;
#print join("\t",@notN)."\n"; 
#			print Dumper(\%na_stat); die;
		}
	}
	close IN;
}

#print Dumper(\%stat); die;

open OUT ,"> $od/$prefix.index.xls";
my @index2 = qw/Ne Ho He nei shi pic/;
print OUT join("\t","sample","Number","Na",@index2,"F_index","maf","PLN","PPL","loci")."\n"; ###PLN => poly loci; PPL => PLN/loci;
for my $g (sort keys %genotype){
	next if ($num{$g} <= 1);
	my $out = "$g\t$stat{$g}{sample_count}";
	my $Na = sprintf("%.4f",$na_stat{$g}{na}/$na_stat{$g}{num});
	$out .= "\t$Na";
	for my $index (@index){
		$out .= "\t".sprintf("%.4f",$stat{$g}{$index}/$stat{$g}{num});
	}
	my $Hoe = $stat{$g}{He} ? $stat{$g}{Ho}/$stat{$g}{He} : 1;
	my $f_index = sprintf("%.4f", 1 - $Hoe);
	$stat{$g}{non_poly} = $stat{$g}{non_poly} ? $stat{$g}{non_poly} : 0;
	my $pln = $stat{$g}{num} - $stat{$g}{non_poly};
	my $ppl = sprintf("%.4f", $pln/$stat{$g}{num});
#	my $Na = sprintf("%.4f",$na_stat{$g}{na}/$na_stat{$g}{num});
	my $maf = $pln > 0 ? sprintf("%.4f",$stat{$g}{maf}/$pln) : "NA";
	$out .= "\t".join("\t",$f_index,$maf,$pln,$ppl,$stat{$g}{num});
	print OUT "$out\n";
}
close OUT;


system("cp $Bin/manual.docx $outdir") if (-s "$Bin/manual.docx");

my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);

###############Subs@##############################
sub deal_geno{
	my ($in,$od,$key,$group,$flag,$genotype) = @_;
	my (%out,%index,%id); 
	if ($in =~ /gz$/){
		open IN, "gzip -dc $in |"; 
	}else{
		open IN, $in or die $!;
	}
	while(<IN>){
		chomp;
		next if (/^##/);
		if ($.==1 || $_=~/^#chr/ || $_=~/^#CHROM/){
			&deal_head($_,$group,$flag,\%index,\%id);
			for my $gg (keys %index){
				open my $out_geno, ">$od/$key.$gg.genotype";
				print $out_geno join("\t", "Chr\tLoci", @{$id{$gg}})."\n";
				$out{$gg} = $out_geno;
				$$genotype{$gg} = "$od/$key.$gg.genotype";
			}
		}else{
			my @x = split/\t/;
			my $loc = join("\t",@x[0,1]);
			for my $gg (keys %index){
				my $out = $out{$gg};
				if ($flag eq "geno"){
					print $out join("\t",$x[0],$x[1],@x[@{$index{$gg}}])."\n";
				}elsif ($flag eq "vcf"){
					if (length($x[3]) == 1 && length($x[4]) == 1){
						my $v2g = &vcf2geno($loc,"$x[3] $x[4]",$x[8],@x[@{$index{$gg}}]);
						print $out $v2g."\n";
					}
				}else{
					if ($x[5] eq "SNP"){
						my $s2g = &GT2geno($loc,"$x[3] $x[4]",@x[@{$index{$gg}}]);
						print $out $s2g."\n";
					}
				}
			}
		}
	}
	close IN;
	close $_ for (values %out);
}

sub GT2geno{
	my ($loc,$base,@lin) =@_;
	my @base=split /\s/,$base;
	for my $gt (@lin){
		$gt=~s/[,\/]/|/;
		if($gt eq "0|0"){
			$loc .= "\t".$base{"$base[0] $base[0]"};
		}elsif($gt eq "0|1"){
			$loc .= "\t".$base{"$base[0] $base[1]"};
		}elsif($gt eq "1|1"){
			$loc .= "\t".$base{"$base[1] $base[1]"};
		}else{
			$loc .= "\tN";
		}
	}
	return $loc;
}

sub vcf2geno{
	my ($loc,$base,$info,@lin) =@_;
	my @info = split/:/,$info;
	my $N = "xxx" ;
	for (my $i=0; $i<= $#info; $i++){
		if ($info[$i] eq "GT"){
			$N = $i; 
			last;
		}
	}
	if ($N eq "xxx"){
		info("ERROR: vcf file has no GT information...");
		exit;
	}
	my @out = map{$_=~s/\//|/;($_ eq ".|.") ? ".|." : (split/:/,$_)[$N];} @lin ;
	return &GT2geno($loc,$base,@out);
}

sub deal_head{
	my ($in,$gg,$flag,$index,$id) = @_;
	my @head = split/\t/,$in;
	my $i = 0;
	if ($flag eq "vcf"){
		@head = @head[9..$#head];
		$i = 9;
	}elsif ($flag eq "snp"){
		my $m; 
		for ($m=8;$m <= $#head; $m++){
			if ($head[$m] =~ /_type$/){
				next;
			}else{
				last;
			}
		}
		@head = @head[8..$m-1];
		$i = 8; 
	}else{
		@head = @head[2..$#head];
		$i = 2;
	}
	my %group = %$gg;
	for my $h (@head){
		$h =~ s/_type$//;
		if (exists $group{all}){
			push @{$$id{all}},$h;
			push @{$$index{all}}, $i;
		}else{
			for my $gg (keys %group){
				for my $sample (@{$group{$gg}}){
					if ($h eq $sample){
						push @{$$id{$gg}},$h;
						push @{$$index{$gg}},$i;
						push @{$$id{all}},$h;
						push @{$$index{all}}, $i;
					}
				}
			}
		}
		$i++;
	}
}

sub deal_group{
	my $group = shift;
	my %g;
	if ($group eq 1){
		$g{all} = 1;
	}else{
		open GG, $group or die $!;
		while(<GG>){
			chomp;
			my ($s,$g) = split/\t/;
			push @{$g{$g}},$s;
			$num{$g} += 1;
			$num{all}+= 1;
		}
		close GG;
	}
	return %g;
}

sub deal_allele{
	my ($arr,$g,$het,$loc) = @_;
	my $allele_line = join("", @$arr); 
	my $number = scalar(@$arr);
	my ($total_allele,$allele_freq_squa_sum,$log_freq_sum,$allele_n,$non_poly) = (0,0,0,0);
	my (%allele,%allele_num,%allele_freq,@key_array,@alleleF);
	map{$allele_num{$_}++; $total_allele ++} split(//, $allele_line); 
	foreach my $key (keys %allele_num){
		$allele_freq{$key} = $allele_num{$key}/$total_allele;
		$allele_freq_squa_sum +=  $allele_freq{$key}*$allele_freq{$key};
		$log_freq_sum +=  $allele_freq{$key}*log($allele_freq{$key});
		push @key_array, $key;
		$allele_n++;
	}
##### MAF stat;
	foreach my $key (sort{$allele_freq{$b} <=> $allele_freq{$a}} keys %allele_freq){
		push(@alleleF, $allele_freq{$key});
	}
	if (scalar @alleleF == 1){
		$stat{$g}{non_poly} ++;   ###非多态性；
		$alleleF[1] = 0;
		$allele{Na}=1;
	}else{
		$allele{Na}=2;
	}
	$stat{$g}{maf} += $alleleF[1];
	$allele{ae} = sprintf("%.4f",1/$allele_freq_squa_sum);    #####
	$allele{Ho} = sprintf("%.4f",$het/scalar @$arr);          ##### Ho : heterozygosis_observed;
	$allele{He} = sprintf("%.4f",1 - $allele_freq_squa_sum);  ##### He : heterozygosis_expected;
	$allele{nei} =sprintf("%.4f",2*(scalar @$arr)*(1 - $allele_freq_squa_sum)/(2*(scalar @$arr) - 1));   ### nei_diversity
	$allele{shi} =sprintf("%.4f",-$log_freq_sum);             ##### Shannon_Wiener_index; 

	my $multi_sum = 0;
	for (my $i = 0; $i <= scalar(keys %allele_num)-2; $i++) {
		for (my $j = $i+1; $j <= $#key_array; $j++) {
			if ($i == 0) {
				$multi_sum = 2*($allele_freq{$key_array[$i]}*$allele_freq{$key_array[$j]})*($allele_freq{$key_array[$i]}*$allele_freq{$key_array[$j]});
			} else {
				$multi_sum = $multi_sum + 2*($allele_freq{$key_array[$i]}*$allele_freq{$key_array[$j]})*($allele_freq{$key_array[$i]}*$allele_freq{$key_array[$j]});
			}
		}
	}
	$multi_sum = 0 if ($allele{ae} == 1);
	$allele{pic} = sprintf("%.4f",$allele{He} - $multi_sum);   ### pic :  Polymorphysm_information_content;
	if (defined $run){
		my $o_loci = $loci{$g};
		print $o_loci join("\t",$loc,$number,$alleleF[1],(map{exists $allele{$_} ? $allele{$_} : 0} @index),$allele{Na})."\n";
	}
	map{$stat{$g}{$_} += (exists $allele{$_}) ? $allele{$_} : 0 } @index;
	$stat{$g}{num}++; 
}
 
###############################################################################################

sub showTime{
	my ($text) = @_;
	my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
	my $format_time = sprintf("[%d-%.2d-%.2d %.2d:%.2d:%.2d]",$year+1900,$mon+1,$mday,$hour,$min,$sec);
	print STDERR "$format_time $text\n";
}

sub info{
	my ($text) = @_;
	&showTime($text);
}


sub sub_format_datetime #Time calculation subroutine
{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Runtime # &Runtime($BEGIN);
{
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total $programe_dir elapsed time : [",&sub_time($t),"]\n";
}
sub sub_time
{
	my ($T)=@_;chomp $T;
	my $s=0;my $m=0;my $h=0;
	if ($T>=3600) {
		my $h=int ($T/3600);
		my $a=$T%3600;
		if ($a>=60) {
			my $m=int($a/60);
			$s=$a%60;
			$T=$h."h\-".$m."m\-".$s."s";
		}else{
			$T=$h."h-"."0m\-".$a."s";
		}
	}else{
		if ($T>=60) {
			my $m=int($T/60);
			$s=$T%60;
			$T=$m."m\-".$s."s";
		}else{
			$T=$T."s";
		}
	}
	return ($T);
}

sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}


sub help
{
	print <<"	Usage End.";
	Description:
	Data    : $Write_Date
	Version : $ver
	function: ......
	Usage:
		[snplist or genotype or vcffile only choose one]
		-i		snplist file		must be given
		-t		genotype file		[table format: chr loc samples ]
		-v		vcf file
		
		-g		group list		[format: sample group; noheader];
		-od		out directory				must be given
		-k		prefix
		-r		output index of each loci 
		
		-s 		genome size;
		-nm		if excu Nm_fst_stat.pl
		-h		Help				document
	
	e.g. perl $0 -i snp.avinput -g group.list -od xxx -k xxx
	     perl $0 -t snp.genotype -g group.list -od xxx -k xxx
	     perl $0 -v snp.vcf -g group.list -od xxx -k xxx
	     perl $0 -i snp.avinput -g group.list -od xxx -k xxx -r   ### output index of each loci
	Usage End.
	exit;
}


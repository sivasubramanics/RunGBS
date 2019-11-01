#!/usr/bin/env perl
# Code: useGBS.pl <option>
# Contact: s.sivasubramani@cgiar.org
# Description: Wrapper to run GBS or ddRAD data.


use Getopt::Long qw(GetOptions);
use Getopt::Long qw(:config no_ignore_case);
use Getopt::Long  qw(:config bundling);

$usage="
	USAGE: $0 --se/--pe -r <reference.fa> -k <GBS_KeyFile.txt> -i <fastq_dir_path>

	Options:
	--se|S	- Single-end data
	--pe|P	- Paired-end data
	--reference|-r	- Reference sequence file. (.fasta or .fa)
	--keyfile|-k	- Path for Keyfile.
	--fastq|-i	- Path for fastq directory (fastq files should named as FLOWCELLIF_LANE.fastq.gz)

	Fastq processing:
	--adapter|-a	- Adapter sequence to use for cutadapt tool. (default: GAGATCGGAAGAGCGGG)
	--minlength|-m	- Minimum read length to be considered. (default: 50)

	Performance:
	--threads|-p	- Number of threads/cores to be used for bowtie2/samtools/bcftools (default: 2)
	--npcore|-n	- Number of threads/cores to be used for parallel (default: 2)
	
	Tools:
	--sabre	- Path for sabre executable (default: \$PATH/sabre)
	--cutadapt	- Path for cutadapt executable (default: \$PATH/cutadapt)
	--bowtie2	- Path for bowtie2 executable (default: \$PATH/bowtie2)
	--samtools	- Path for samtools executable (samtools-1.6 or higher) (default: \$PATH/samtools)
	--bcftools	- Path for bcftools executable (bcftools-1.6 or higher) (default: \$PATH/bcftools)
	--bamtools	- Path for bcftools executable (default: \$PATH/bamtools)
";


$help=0;
$NPCORE=2;
$NTH=2;
$sabre="sabre";
$cutadapt="cutadapt";
$bowtie2="bowtie2";
$samtools="samtools";
$bcftools="bcftools";
$bamtools="bamtools";
$NTH=2;
$MINLN=50;
@lanes=();
$ADAPTER="GAGATCGGAAGAGCGGG";
$SE=0;
$PE=0;
%samples = ();
$flag = 0;


GetOptions(
	"reference|r=s" => \$REFERENCE,
	"keyfile|k=s" => \$KeyFile,
	"fastq|i=s" => \$FASTQ,
	"adapter|a:s" => \$ADAPTER,
	"minlength|m:s" => \$MINLN,
	"threads|p:i" => \$NTH,
	"npcore|n:i" => \$NPCORE,
	"help|h" => \$help,
	"bowtie2:s" => \$bowtie2,
	"samtools:s" => \$samtools,
	"bcftools:s" => \$bcftools,
	"cutadapt:s" => \$cutadapt,
	"sabre:s" => \$sabre,
	"se|S" => \$SE,
	"pe|P" => \$PE
	)
or die "$usage\n";

if($help==1 || $REFERENCE eq ""){
	die "$usage\n";
}

if($SE eq $PE){
	die "Invalid data set option. Please mention either single-end data or paired-end data. $usage\n";
}

if($REFERENCE=~/(.*).fasta$/ || $REFERENCE=~/(.*).fa$/){
	$REF_BASE=$1;
}

chomp($ALLCORE=`cat /proc/cpuinfo | grep processor | tail -1 | awk \'{print \$NF}\'`);
$KeyFileSorted=$KeyFile.".sorted";
$DATADIR=$FASTQ."/Analysis";

if(! -d $DATADIR){
	print "Creating Analyis  directory... ($DATADIR)\n";
	$exit_status = mkdir($DATADIR);
	checkLog($exit_status,"makeDataDirectory");
}

if(-f "$DATADIR/check_point.txt"){
	open(CHECK, "$DATADIR/check_point.txt");
	while (<CHECK>) {
		chomp();
		$check_point{$_}=1;
	}
	close CHECK;
}

print "Creating check point file... ($DATADIR/check_point.txt)\n";
open(CHECK, ">>$DATADIR/check_point.txt");

$sortCMD = "sort -k1,1 -k2,2 $KeyFile > $KeyFileSorted";
print "Sorting Keyfile for FlowCell and Lane... ($KeyFileSorted)\n";
$exit_status = system($sortCMD);
checkLog($exit_status,"sortKeyFile");


open(KEY, $KeyFileSorted);
while(<KEY>){
	chomp();
	$line=$_;
	@arr = split("\t", $line);
	if($arr[0] eq "Flowcell"){
		next;
	}
	$out=$arr[0]."_".$arr[1];
	if($pout ne $out){
		push(@lanes, $out);
		$out_file=$FASTQ."/".$out.".barcodes";
		print "Creating Barcode file fr sabre... ($out_file)\n";
		open($OF, ">$out_file");
	}
	$fq{$arr[3]}{$arr[0]}{$arr[1]}=1;
	$samples{$arr[3]}++;
	if($SE==1){
		print $OF $arr[2],"\t",$FASTQ,"/",$arr[0],"_",$arr[1],"/",$arr[3],"_n",$samples{$arr[3]},".fq\n";
	}
	if($PE==1){
		print $OF $arr[2],"\t",$FASTQ,"/",$arr[0],"_",$arr[1],"/",$arr[3],"_n",$samples{$arr[3]},"_1.fq","\t",$FASTQ,"/",$arr[0],"_",$arr[1],"/",$arr[3],"_n",$samples{$arr[3]},"_2.fq\n";
	}
	$pout = $out;
}
close KEY;


if(! exists $check_point{"SABRE"}){
	foreach $lane(@lanes){
		$barcodeFile = $FASTQ."/".$lane.".barcodes";
		if($SE==1){
			$fastqFile = $FASTQ."/".$lane."_fastq.gz";
			if(! -f $fastqFile && $fastqFile ne ""){
				print "WARNING: $fastqFile file doesn't exist. Skipping\n";
			}
			else{
				mkdir("$FASTQ/$lane");
				$sabreCMD="\"$sabre se -f $fastqFile -b $barcodeFile -u $FASTQ/$lane/unknown.fq \$\> $FASTQ/$lane/sabre.log\"";
				push(@sabreCMDs, $sabreCMD);
			}
		}
		if($PE==1){
			$fastqFile1 = $FASTQ."/".$lane."_1.fastq.gz";
			$fastqFile2 = $FASTQ."/".$lane."_2.fastq.gz";
			if((! -f $fastqFile1) || (!-f $fastqFile2)){
				print "WARNING: $fastqFile1 or $fastqFile1 file doesn't exist. Skipping\n";
			}
			else{
				mkdir("$FASTQ/$lane");
				$sabreCMD="\"$sabre pe -f $fastqFile1 -r $fastqFile2 -b $barcodeFile -u $FASTQ/$lane/unknown_1.fq -w $FASTQ/$lane/unknown_2.fq \$\> $FASTQ/$lane/sabre.log\"";
				push(@sabreCMDs, $sabreCMD);
			}
		}
	}
	$sabreCMD="parallel -j $NPCORE \:\:\: ".join(" ", @sabreCMDs);
	print $sabreCMDs,"\n";
	print "Processing Sabre...\n";
	$exit_status = system($sabreCMD);
	checkLog($exit_status,"sabre");
	print CHECK "SABRE\n";
}

if(! exists $check_point{"CONCATENATE"}){
	print "Concatenating fastq files...\n";
	foreach $ID(sort keys %samples){
		if($SE==1){
			$file="";
			for $lane(@lanes){
				for ($i=1;$i<=$samples{$ID};$i++){
					$fq = $FASTQ."/".$lane."/".$ID."_n".$i.".fq";
					$fqgz = $FASTQ."/".$lane."/".$ID."_n".$i.".fq.gz";
					if(-f $fq){
						$file.= " ".$fq;
					}
					elsif(-f $fqgz){
						$ungzipCMD = "pigz -p $ALLCORE -d $fqgz";
						$exit_status = system($ungzipCMD);
						checkLog($exit_status,"decompressing");
						$file.= " ".$fq;	
					}
				}
			}
			$fqfile=$DATADIR."/".$ID.".fq.gz";
			if($file ne ""){
				$catCMD = "cat $file | pigz --fast -p $ALLCORE >> $fqfile";
				$exit_status = system($catCMD);
				checkLog($exit_status,"concatenating");
				$gzipCMD = "pigz --fast -p $ALLCORE $file";
				$exit_status = system($gzipCMD);
				checkLog($exit_status,"compressing");
			}
		}
		if($PE==1){
			$file1="";
			$file2="";
			for $lane(@lanes){
				for ($i=1;$i<=$samples{$ID};$i++){
					$fq1 = $FASTQ."/".$lane."/".$ID."_n".$i."_1.fq";
					$fq2 = $FASTQ."/".$lane."/".$ID."_n".$i."_2.fq";
					$fq1gz = $FASTQ."/".$lane."/".$ID."_n".$i."_1.fq.gz";
					$fq2gz = $FASTQ."/".$lane."/".$ID."_n".$i."_2.fq.gz";
					if(-f $fq1){
						$file1.= " ".$fq1;
					}
					elsif(-f $fq1gz){
						$ungzipCMD = "pigz -p $ALLCORE -d $fq1gz";
						$exit_status = system($ungzipCMD);
						checkLog($exit_status,"decompressing");
						$file.= " ".$fq1;
					}
					if(-f $fq2){
						$file2.= " ".$fq2;
					}
					elsif(-f $fq2gz){
						$ungzipCMD = "pigz -p $ALLCORE -d $fq2gz";
						$exit_status = system($ungzipCMD);
						checkLog($exit_status,"decompressing");
						$file.= " ".$fq2;
					}
				}
			}
			$fqfile1=$DATADIR."/".$ID."_1.fq.gz";
			$fqfile2=$DATADIR."/".$ID."_2.fq.gz";
			if($file1 ne "" && $file2 ne ""){
				$catCMD1 = "cat $file1 | pigz --fast -p $ALLCORE >> $fqfile1";
				$exit_status = system($catCMD1);
				checkLog($exit_status,"concatenating");
				$gzipCMD = "pigz --fast -p $ALLCORE $file1";
				$exit_status = system($gzipCMD);
				checkLog($exit_status,"compressing");

				$catCMD2 = "cat $file2 | pigz --fast -p $ALLCORE >> $fqfile2";
				$exit_status = system($catCMD2);
				checkLog($exit_status,"concatenating");
				$gzipCMD = "pigz --fast -p $ALLCORE $file2";
				$exit_status = system($gzipCMD);
				checkLog($exit_status,"compressing");
			}	
		}
	}
	print CHECK "CONCATENATE\n";
}

if(! exists $check_point{"CUTADAPT"}){
	if($SE==1){
		$cutadaptCMD="parallel -j $NPCORE $cutadapt --quiet -a $ADAPTER -m $MINLN -o {}.fastq.gz {}.fq.gz \'2\>\' {}.cutadapt.log \:\:\: \$(ls -1 $DATADIR/*.fq.gz | sed 's/.fq.gz//')";
	}
	if($PE==1){
		$cutadaptCMD="parallel -j $NPCORE $cutadapt --quiet -a $ADAPTER -A $ADAPTER -m $MINLN -o {}_1.fastq.gz -p {}_2.fastq.gz {}_1.fq.gz {}_2.fq.gz \'2\>\' {}.cutadapt.log \:\:\: \$(ls -1 $DATADIR/*_1.fq.gz | sed 's/_1.fq.gz//')";
	}
	print "Processing Cutadapt...\n";
	print $cutadapt,"\n";
	$exit_status = system($cutadaptCMD);
	checkLog($exit_status,"cutadapt");
	print CHECK "CUTADAPT\n";
}

if(! -f $REF_BASE.".1.bt2"){
	print "Bowtie2 index does not exist. Creating one.\n";
	$bowtie2idxCMD="$bowtie2\-build $REFERENCE $REF_BASE";
	print "Creating reference bowtie2 index...\n";
	$exit_status = system($bowtie2idxCMD);
	checkLog($exit_status,"bowtie2Indexing");
}

if(! -f $REFERENCE.".fai"){
	print "fasta index does not exist. Creating one.\n";
	$fastaidxCMD="$bowtie2\-build $REFERENCE $REF_BASE";
	print "Creating reference fasta index...\n";
	$exit_status = system($fastaidxCMD);
	checkLog($exit_status,"fastaIndexing");
}

open(REFIDX, $REFERENCE.".fai");
while (<REFIDX>) {
	chomp();
	@arr = split("\t", $_);
	$chrFileName = $REF_BASE."_".$arr[0].".fa";
	push(@refereceChrs, $arr[0]);
	if(! -f $chrFileName){
		$flag = 1;
	}
}
close REFIDX;

if($flag == 1){
	open(REF, $REFERENCE);
	while (<REF>) {
		chomp();
		if($_=~/^>(\S+)/){
			$chrName = $1;
			$chrFileName = $REF_BASE."_".$chrName.".fa";
			close CHR;
			open(CHR, ">$chrFileName");
			print CHR $_,"\n";
		}
		else{
			print CHR $_,"\n";	
		}
	}
	close CHR;
	close REF;
}

if(! exists $check_point{"BOWTIE2"}){
	open(BOWTIESHELL, ">$DATADIR/runBowtie2.sh");
	if($SE==1){

		foreach $sample(keys %samples){
			$tempFq = $DATADIR."/".$sample.".fastq.gz";
			if(-f $tempFq){
				print BOWTIESHELL "$bowtie2 --end-to-end -N 1 -p $NTH --rg-id $sample --rg PL:ILLUMINA --rg SM:$sample -x $REF_BASE -U $DATADIR/$sample.fastq.gz 2> $DATADIR/$sample.log | $samtools view -T $REFERENCE -hu - | $samtools sort --reference $REFERENCE -o $DATADIR/$sample.sort.bam 2> $DATADIR/$sample.log\n";	
			}
			
		}

		# $bowtie2CMD="parallel -j $NPCORE $bowtie2 --end-to-end -p $NTH --rg-id {/.} --rg \"PL:ILLUMINA\" --rg \"SM:{/.}\" -x $REF_BASE -U {}.fastq -S {}.sam \'2\>\' {}.log \:\:\: \$(ls -1 $DATADIR/*.fastq | sed 's/.fastq//')";
	}
	if($PE==1){
		# $bowtie2CMD="parallel -j $NPCORE $bowtie2 --end-to-end -X 500 -p $NTH --rg-id {/.} --rg \"PL:ILLUMINA\" --rg \"SM:{/.}\" -x $REF_BASE -1 {}_1.fastq -2 {}_2.fastq -S {}.sam \'2\>\' {}.log \:\:\: \$(ls -1 $DATADIR/*_1.fastq | sed 's/_1.fastq//')";
		foreach $sample(keys %samples){
			$tempFq1 = $DATADIR."/".$sample."_1.fastq.gz";
			$tempFq2 = $DATADIR."/".$sample."_2.fastq.gz";
			if(-f $tempFq1 && -f $tempFq2){
				print BOWTIESHELL "$bowtie2 --end-to-end -N 1 -p $NTH --rg-id $sample --rg PL:ILLUMINA --rg SM:$sample -x $REF_BASE -1 $DATADIR/$sample_1.fastq.gz -2 $DATADIR/$sample_2.fastq.gz 2> $DATADIR/$sample.log | $samtools view -T $REFERENCE -hu - | $samtools sort --reference $REFERENCE -o $DATADIR/$sample.sort.bam 2> $DATADIR/$sample.log";
			}
		}

		# $bowtie2CMD="parallel -j $NPCORE $bowtie2 --end-to-end -N 1 -p $NTH --rg-id {/.} --rg \"PL:ILLUMINA\" --rg \"SM:{/.}\" -x $REF_BASE -1 {}_1.fastq -2 {}_2.fastq | $samtools view -T $REFERENCE -hu - | samtools sort --reference $REFERENCE -o {}.sort.bam \'2\>\' {}.log \:\:\: \$(ls -1 $DATADIR/*.fastq | sed 's/.fastq//')";
	}
	close BOWTIESHELL;
	$bowtie2CMD = "parallel -j $NPCORE < $DATADIR/runBowtie2.sh";
	print "Bowtie2 mapping is in progress...\n";
	print "$bowtie2CMD\n";
	$exit_status = system($bowtie2CMD);
	checkLog($exit_status, "bowtie2");
	print CHECK "BOWTIE2\n";
}

=head
if(! exists $check_point{"SAM2BAM"}){
	if($SE==1){
		$sam2bamCMD="parallel -j $NPCORE $samtools view -T $REFERENCE -S {}.sam -b -o {}.bam \:\:\: \$(ls -1 $DATADIR/*.fastq | sed 's/.fastq//')";
	}
	if($PE==1){
		$sam2bamCMD="parallel -j $NPCORE $samtools view -T $REFERENCE -S {}.sam -b -o {}.bam \:\:\: \$(ls -1 $DATADIR/*_1.fastq | sed 's/_1.fastq//')";
	}
	print "SAM to BAM is in progress...\n";
	$exit_status = system($sam2bamCMD);
	checkLog($exit_status, "sam2bam");
	print CHECK "SAM2BAM\n";
}

if(! exists $check_point{"SORTBAM"}){
	if($SE==1){
		$sortbamCMD="parallel -j $NPCORE $samtools sort --threads $NTH {}.bam -o {}.sort.bam \:\:\: \$(ls -1 $DATADIR/*.fastq | sed 's/.fastq//')";
	}
	if($PE==1){
		$sortbamCMD="parallel -j $NPCORE $samtools sort --threads $NTH {}.bam -o {}.sort.bam \:\:\: \$(ls -1 $DATADIR/*_1.fastq | sed 's/_1.fastq//')";
	}
	print "BAM sorting is in progress...\n";
	$exit_status =  system($sortbamCMD);
	checkLog($exit_status, "sortBAM");
	print CHECK "SORTBAM\n";
}
=cut

if(! exists $check_point{"SORTBAMIDX"}){
	if($SE==1){
		$sortbamidxCMD="parallel -j $ALLCORE $samtools index {}.sort.bam \:\:\: \$(ls -1 $DATADIR/*.fastq.gz | sed 's/.fastq.gz//')";
	}
	if($PE==1){
		$sortbamidxCMD="parallel -j $ALLCORE $samtools index {}.sort.bam \:\:\: \$(ls -1 $DATADIR/*_1.fastq.gz | sed 's/_1.fastq.gz//')";
	}
	print "BAM indexing is in progress...\n";
	$exit_status =  system($sortbamidxCMD);
	checkLog($exit_status, "sortBAMindex");
	print CHECK "SORTBAMIDX\n";
}

if(! exists $check_point{"SPLITBAM"}){
	open(BAMSPLIT, ">$DATADIR/runSplitBam.sh");
	foreach $sample(keys %samples){
		$tempBam = $DATADIR."/".$sample.".sort.bam";
		if(-f $tempBam){
			print BAMSPLIT "$bamtools split -refPrefix \'\' -in $tempBam -reference\n";
		}
	}
	close BAMSPLIT;
	$splitBamCMD = "parallel -j $ALLCORE < $DATADIR/runSplitBam.sh";
	print "BAM splitting is in progress...\n";
	print "$splitBamCMD\n";
	$exit_status =  system($splitBamCMD);
	checkLog($exit_status, "splitBAM");
	print CHECK "SPLITBAM\n";
}

if(! exists $check_point{"SPLITBAMIDX"}){
	open(BAMSPLITIDX, ">$DATADIR/runSplitBamIndex.sh");
	foreach $sample(keys %samples){
		$tempBam = $DATADIR."/".$sample.".sort.bam";
		if(-f $tempBam){
			print BAMSPLITIDX "$samtools index $tempBam\n";
		}
	}
	close BAMSPLITIDX;

	$splitBamIndexCMD = "parallel -j $ALLCORE < $DATADIR/runSplitBam.sh";
	print "BAM indexing is in progress...\n";
	print "$splitBamIndexCMD\n";
	$exit_status =  system($splitBamIndexCMD);
	checkLog($exit_status, "splitBAMIndex");
	print CHECK "SPLITBAMIDX\n";
}

if(! exists $check_point{"SORTBAMLIST"}){
	foreach $chr(@refereceChrs){
		$bamList="ls -c1 $DATADIR/*.sort.$chr.bam > $DATADIR/sortBamList_$chr";	
		$exit_status =  system($bamList);
		checkLog($exit_status, "makeSortBamList_$chr");
	}
	# $bamList="ls -c1 $DATADIR/*.sort.bam > $DATADIR/sortBamList";
	# print "Creating Bam list file ($DATADIR/sortBamList)...\n";
	# $exit_status =  system($bamList);
	# checkLog($exit_status, "makeSortBamList");
	print CHECK "SORTBAMLIST\n";
}

# $rmdupCMD="parallel -j $NPCORE $samtools rmdup -s {}.sort.bam {}.rmdup.bam \:\:\: \$(ls -1 $DATADIR/*.fastq | sed 's/.fastq//')";
# print $rmdupCMD,"\n";
# $exit_status =  system($rmdupCMD);
# checkLog($exit_status, "rmdup");

# $rmdupbamidxCMD="parallel -j $NPCORE $samtools index {}.rmdup.bam \:\:\: \$(ls -1 $DATADIR/*.fastq | sed 's/.fastq//')";
# print $rmdupbamidxCMD,"\n";
# $exit_status =  system($rmdupbamidxCMD);
# checkLog($exit_status, "rmdupBAMindex");

# $bamList="ls -c1 $DATADIR/*.rmdup.bam > $DATADIR/rmdupBamList";
# print $bamList,"\n";
# $exit_status =  system($bamList);
# checkLog($exit_status, "makeRMdupBamList");

# $varcallCMD="$bcftools mpileup -a DP,AD -O u --threads $NTH -f $REFERENCE --bam-list $DATADIR/rmdupBamList | bcftools call --threads $NTH -mv -O z -o $DATADIR/var.raw.rmdup.vcf.gz";
# print $varcallCMD,"\n";
# $exit_status =  system($varcallCMD);
# checkLog($exit_status, "VariantCallingRMduptBAM");

if(! exists $check_point{"VARCALLING"}){
	open(VARCALL, ">$DATADIR/runVarCalling.sh");
	foreach $chr(@refereceChrs){
		$chrFileName = $REF_BASE."_".$chr.".fa";
		$varcallCMD="$bcftools mpileup -a DP,AD -O u --threads $NTH -f $chrFileName --bam-list $DATADIR/sortBamList_$chr | $bcftools call --threads $NTH -mv -O z -o $DATADIR/var.$chr.raw.sort.vcf.gz";
		print VARCALL $varcallCMD,"\n";
	}
	close VARCALL;
	$varcallCMD = "parallel -j $NPCORE < $DATADIR/runVarCalling.sh";
	print "Variant calling is in progress...\n";
	print "$varcallCMD\n";
	$exit_status =  system($varcallCMD);
	checkLog($exit_status, "VariantCallingSortBAM");
	print CHECK "VARCALLING\n";

	# $varcallCMD="$bcftools mpileup -a DP,AD -O u --threads $ALLCORE -f $REFERENCE --bam-list $DATADIR/sortBamList | $bcftools call --threads $NTH -mv -O z -o $DATADIR/var.raw.sort.vcf.gz";
	# print "Variant calling is in progress ($DATADIR/var.raw.sort.vcf.gz)...\n";
	# print $varcallCMD,"\n";
	# $exit_status =  system($varcallCMD);
	# checkLog($exit_status, "VariantCallingSortBAM");
	# print CHECK "VARCALLING\n";
}

if(! exists $check_point{"VARFILTERING"}){
	open(VARFILTER, ">$DATADIR/runVarFiltering.sh");
	foreach $chr(@refereceChrs){
		$varFilterCMD="$bcftools view --types snps -O z -o $DATADIR/var.$chr.raw.sort.snp.vcf.gz $DATADIR/var.$chr.raw.sort.vcf.gz";
		
		print VARFILTER $varFilterCMD,"\n";
	}
	close VARFILTER;
	$varcallCMD = "parallel -j $NPCORE < $DATADIR/runVarFiltering.sh";
	print "Variant Filtering is in progress...\n";
	$exit_status =  system($varFilterCMD);
	checkLog($exit_status, "VariantFilteringForSNPs");
	print CHECK "VARFILTERING\n";
}

close CHECK;

sub checkLog(){
	$exit_sts = shift @_;
	$process = shift @_;
	if($exit_sts > 1){
		die "Something wrong with \"$process\" step. Quiting.\n";
	}
}

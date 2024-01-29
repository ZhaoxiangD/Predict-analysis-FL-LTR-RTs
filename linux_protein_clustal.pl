

open ltr_arm_fasta_line,"$ARGV[0]" || die "Input file ltr_arm_fasta_line cannot be opened.\n"; 

foreach $a (<ltr_arm_fasta_line>) {
	chomp($a);
	if ($a=~/^>(\S+)\t+(\S+)/) {
		$i++;
		$seq_id[$i]=$1;
		$seq{$seq_id[$i]}=$2;
		next;
	}
	if ($a=~/^>(\S+)\s+(\S+)/) {
		$i++;
		$seq_id[$i]=$1;
		$seq{$seq_id[$i]}=$2;
		next;
	}
}
$i_max=$i;
print STDERR ("i_max= $i_max\n");

$file=$ARGV[0].".clustalw";
system "rm $file";

for ($i=1; $i<=$i_max; $i=$i+2) {
	$clustalw_input_file=$file.".current_seq";
	open clustalw_input, ">$clustalw_input_file" || die "clustalw input file cannot be open.\n";
	print clustalw_input (">$seq_id[$i]\n$seq{$seq_id[$i]}\n>$seq_id[$i+1]\n$seq{$seq_id[$i+1]}\n");
	close clustalw_input;
	system "/home/software/clustalw-2.0.12/src/clustalw2 -INFILE=$clustalw_input_file -TYPE=PROTEIN -OUTORDER=INPUT";
	
	$clustalw_output_file=$file.".aln";
	open clustalw_output, "$clustalw_output_file" || die "Input file clustalw_output cannot be opened.\n";
	@text1=<clustalw_output>;
	close clustalw_output;
	
	open output, ">>$ARGV[0].clustalw" || die "Output file cannot be opened.\n";
	foreach $a (@text1) {
		if ($a=~/(\S+)\s+(\S+)$/) {
			chomp($a);
			@tmp=split(" ",$a);
			$seq_ltr_arm{$tmp[0]} = $seq_ltr_arm{$tmp[0]}.$tmp[1];
		}
	}
	print output (">$seq_id[$i]\n$seq_ltr_arm{$seq_id[$i]}\n>$seq_id[$i+1]\n$seq_ltr_arm{$seq_id[$i+1]}\n____________end of seq pair____________\n");
	close output;
}
	
	







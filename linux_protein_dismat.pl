
open ltr_arm_fasta_line_clustalw,"$ARGV[0]" || die "Input file ltr_arm_fasta_line_clustalw cannot be opened.\n"; 

$seq="";
$file=$ARGV[0].".distmat";
system "rm $file";

foreach $a (<ltr_arm_fasta_line_clustalw>) {
	chomp($a);
	if($a eq "") {
		next;
	}
	if ($a =~/____________end of seq pair____________/) {
		$input_file=$file.".distmat_input";
		open distmat_input, ">$input_file" || die "Input file distmat_input file cannot be opened.\n"; 
		@tmp=split(" ",$seq);
		$seq1=$tmp[0];
		$seq{$seq1}=$tmp[1];
		$seq2=$tmp[2];
		$seq{$seq2}=$tmp[3];
		print distmat_input ("$seq1\n$seq{$seq1}\n$seq2\n$seq{$seq2}\n");
		close distmat_input;
		$seq="";
		
		$output_file=$file.".distmat_output";
		system "/home/software/EMBOSS-6.5.7/emboss/distmat -sequence $input_file -protmethod 2 -outfile $output_file";
		
		open distmat_output, "$output_file" || die "Input file distmat_output file cannot be opened.\n"; 
		@text=<distmat_output>;
		close distmat_output;
		
		foreach $b (@text) {
			$b=~s/\s/\t/g;
			$seq_label=substr($seq1,1);
			if ($b=~/\t+(\S+)\t+(\S+)\t+(\S+)\t+(\S+)/ && $b=~/$seq_label/) {
				@tmp=split(/\t+/,$b);
				$ltr_divergence=$tmp[2];
				open output, ">>$ARGV[0].distmat" || die "Input file output file cannot be opened.\n"; 
				print output ("$seq1\t$seq2\t$ltr_divergence\n$seq{$seq1}\n$seq{$seq2}\n\n");
				close output;
				last;
			}
		}
	}
	else {
		$seq=$seq." ".$a;
	}
		
}



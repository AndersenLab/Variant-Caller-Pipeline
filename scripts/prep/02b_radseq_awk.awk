function max(arr)
{
	cur_max = 0;
	for(i in arr) {
		if (arr[i] > cur_max) {
			cur_max = i
		}
	}
	return cur_max
}

function min(arr)
{
	cur_min = 10000;
	for(i in arr) {
		if (arr[i] < cur_min) {
			cur_min = i
		}
	}
	return cur_min
}

function allele_code(allele,ref,alt) {
	if (allele == ref) {
		return "0/0";
	} else if (allele == alt) {
		return "1/1";
	} else {
		return "./.";
	}
}

{ ORS=""; 
	for(item in freq) {
		delete freq[item];
	}
	for ( i=4; i<=NF; i++ ) {
	freq[$i]++;
	}
	ref_allele=max(freq);
	alt_allele=min(freq);
	# Print data
	print "chr" gensub("_.*$","","g",$1); #Chromosome
	print "\t" $2 "\t.\t";
	print ref_allele "\t"; # Ref Allele
	print alt_allele "\t.\tPASS\t.\tGT\t"; # Alt Allele, QUAL, FILTER, INFO, FORMAT
	for ( i=4; i<=NF; i++ ) print allele_code($i,ref_allele, alt_allele) "\t"; print "\n"
}

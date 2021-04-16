#!/bin/bash
#USAGE:
#	crisprviz.sh   [-f fasta filename -optional]
#				[-c clean (BOOLEAN)-optional]
#				[-p run in parallel (BOOLEAN) -optional]
#				[-x split loci (BOOLEAN) -optional]
#				[-r min # of repeats -optional - default = 4]
#				[-b min length of CRIPSR repeats -optional - default = 23]
#				[-e max length of CRISPR repeats -optional - default = 47]
#				[-m min length of CRISPR spacers -optional - default = 26]
#				[-n max length of CRISPR spacers -optional - default 50]
#				[-i include organisms with no detectable CRISPR arrays]
#				[-h show help menu]
#	crisprviz.sh run with no options will process all .fasta files in the directory

while getopts f:cr:b:e:m:n:hpstxi name
do
	case $name in
		p) popt=1;;
		f) fopt=$OPTARG;;
		c) copt=1;;
		r) ropt=$OPTARG;;
		b) bopt=$OPTARG;;
		e) eopt=$OPTARG;;
		m) mopt=$OPTARG;;
		n) nopt=$OPTARG;;
		s) sopt=1;;
		t) topt=1;;
		x) xopt=1;;
		i) iopt=1;;
		*) echo "Invalid arg";;
	esac
done

SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MINCED=${SOURCE_DIR}/java/minced/minced

echo "Beginning crispr locus analysis..."

##setup basic run params
MIN_REPEATS="-minNR 4"
if [ ! -z $ropt ]; then
	MIN_REPEATS="-minNR $ropt"
	echo "Setting min # of repeats: $MIN_REPEATS"
fi

MIN_REPEAT_LEN=""
if [ ! -z $bopt ]; then
	MIN_REPEAT_LEN=" -minRL $bopt"
	echo "Setting min length of repeats: $MIN_REPEAT_LEN"
fi

MAX_REPEAT_LEN=""
if [ ! -z $eopt ]; then
	MAX_REPEAT_LEN=" -maxRL $eopt"
	echo "Setting max length of repeats: $MAX_REPEAT_LEN"
fi

MIN_SPACER_LEN=""
if [ ! -z $mopt ]; then
	MIN_SPACER_LEN=" -minSL $mopt"
	echo "Setting min length of spacers: $MIN_SPACER_LEN"
fi

MAX_SPACER_LEN=""
if [ ! -z $nopt ]; then
	MAX_SPACER_LEN=" -maxSL $nopt"
	echo "Setting max length of spacers: $MAX_SPACER_LEN"
fi


extract_repeats () {
	#check for single file run
	if [ ! -z ${1} ]; then
		gen_repeat_file "${1}"
	else
		CRISPR_FILES=()
		while IFS=  read -r -d $'\0'; do
	    	CRISPR_FILES+=("$REPLY")
		done < <(find . -maxdepth 1 -name '*.crisprs' -print0)
		echo "Total .crisprs files available in directory: ${#CRISPR_FILES[@]}"

		for i in "${CRISPR_FILES[@]}"
		do
			gen_repeat_file "${i}"
		done
	fi

	#remove working tmp files
	rm ${INPUT_FILE}.crisprs_*.tmp
}

gen_repeat_file () {
	grep "^[0-9]" "${1}" | awk '{print $2}' > "${1}_repeats.base.tmp"
	repeatFileName="${1//.crisprs}"
	declare -a repeats_base
	i=1
	while IFS=$'\n' read -r line_data; do
    	repeats_base[i]="${line_data}"
    	((i++))
	done < "${1}_repeats.base.tmp"
	# echo "Repeats base: ${repeats_base[*]}"

	# echo "Generating repeat file for : ${1}"
	# rm "${repeatFileName}_repeats.fa"f
	repeatIterator=1
	repeatIdx=1
	repeatLocus=1
	grep "^[0-9]" "${1}" > "${1}_repeats.tmp"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		# lineGap=$(awk '{print $3}' $line)
		IFS='[ ' read -r -a repeat_lines <<< "$line"

		# echo "writing to file: ${repeats_base[${repeatIdx}]}"
  		echo ">repeat_locus${repeatLocus}_${repeatIdx}" >> "${repeatFileName}_repeats.fa"
  		echo "${repeats_base[${repeatIterator}]}" >> "${repeatFileName}_repeats.fa"

  		# echo "${line} ${#repeat_lines[@]}"
  		# echo "${#repeat_lines[@]} ${repeatIdx}"
    	if [ "${#repeat_lines[@]}" -lt 2 ]; then #write to file and start over
    		((repeatLocus++))
    		repeatIdx=0
    	fi
    	# echo "repeatIdx: $repeatIdx"
    	((repeatIdx++))
    	((repeatIterator++))
	done < "${1}_repeats.tmp"
}

get_spacer_files () {
	SPACER_FILES=()
	while IFS=  read -r -d $'\0'; do
	    SPACER_FILES+=("$REPLY")
	done < <(find . -maxdepth 1 -name '*_spacers.fa' -print0)
	echo "Total spacer files available in directory: ${#SPACER_FILES[@]}"

}

get_repeat_files () {
	REPEAT_FILES=()
	while IFS=  read -r -d $'\0'; do
    	REPEAT_FILES+=("$REPLY")
	done < <(find . -maxdepth 1 -name '*_repeats.fa' -print0)
	echo "Total repeat files available in directory: ${#REPEAT_FILES[@]}"

}

remove_null_results () {
	if [[ -z ${iopt} ]]; then
		echo "Removing null results..."
		find . -maxdepth 1 -empty -delete
	fi
}

#check for single file parameter
INPUT_FILE=""
if [ ! -z $fopt ]; then
	echo "Setting input file: "${fopt}
	INPUT_FILE="$fopt"


	${MINCED} -spacers ${MIN_REPEATS}${MIN_REPEAT_LEN}${MAX_REPEAT_LEN}${MIN_SPACER_LEN}${MAX_SPACER_LEN} ${INPUT_FILE} ${INPUT_FILE}.crisprs
	extract_repeats ${INPUT_FILE}.crisprs
	remove_null_results
fi

echo "Minced parsing completed @: " $(date)

exit 0



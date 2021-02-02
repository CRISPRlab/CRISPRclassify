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
		h) hopt=1;;
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
# CRISPRVIS_ENGINE=${SOURCE_DIR}/crisprviz_engine.py



if [[ -n $hopt && $hopt -gt 0 ]]; then
	echo '''
CRISPRviz <>||<>||<>, a tool for spacer|repeat identification and exploration in full genomes

Usage: crisprviz.sh [OPTIONS]

Options:	-t test installation. Will run the pipeline on a set of test data. Go to localhost:4444 in your browser to see results	
		-p 	parallel processing. Concurrently processes all genome.fasta files in current directory | (recommended for < 50 - 100 files) 	
		-s 	run crisprviz_engine only. Use when spacers|repeats have already been extracted.
		-x 	split loci. Each loci gets listed as a separate row in web tool results
		-f 	.fasta genome filename - file must be in ".fasta" format. *NOTE: Locates file in current directory ONLY | (optional) 
		-c 	clean. remove *.crisprs, *_spacers.fa, *.crisprs.fa prior to execution | (optional)
		-r 	min # of repeats | Default = 4 | (optional)
		-b 	min length of CRIPSR repeats in array | Default = 23 (optional)
		-e 	max length of CRISPR repeats in array | Default = 47 (optional)
		-m 	min length of CRISPR spacers in array | Default = 26 (optional)
		-n 	max length of CRISPR spacers in array | Default = 50 (optional)
		-i 	include organisms with no detecatble CRISPR arrays (Off by default)
		-h 	show this help menu

Examples:	# Standard run: processes all genome.fasta files in current directory in parallel (concurrently), cleans, and splits loci.
			crisprviz.sh -pxc

		# Run for a single genome file, with a min # of repeats of 3, and clean directory
			crisprviz.sh -cx -f genome.fasta -r 3

		# Run for all current directory files with min and max spacer length of 28 and 53 respectively
			crisprviz.sh -px -l 28 -x 53
		'''
	exit 1
fi



# execute_spacer_analysis () {
# 	RUN_FILE=""
# 	if [ ! -z ${1} ]; then
# 		RUN_FILE="-f ${1}"
# 	fi

# 	${CRISPRVIS_ENGINE} ${RUN_FILE} "-p${SOURCE_DIR}"
# }

## run on test data
# if [[ -n $topt && $topt -gt 0 ]]; then
# 	echo "Setting current dir: ${SOURCE_DIR}/../test"
# 	cd ${SOURCE_DIR}/../test
# 	copt=1
# 	fopt=""
# 	sopt=0
# 	popt=1
# 	xopt=1
# fi

## only run spacer compression engine
# if [[ -n $sopt && $sopt -gt 0 ]]; then
# 	SPACER_INPUT_FILE=""
# 	if [ ! -z $fopt ]; then 
# 		SPACER_INPUT_FILE=${fopt}
# 		# must be *_spacers.fa format
# 		execute_spacer_analysis ${SPACER_INPUT_FILE}
# 		exit 1
# 	fi

# 	${CRISPRVIS_ENGINE} "-p${SOURCE_DIR}" 
# 	exit 1
# fi


echo "Beginning crispr locus analysis..."

## clean out existing files from previous runs ##
# if [[ -n $copt && $copt -gt 0 ]]; then
# 	rm ./data/*.crisprs
# 	rm ./data/*_spacers.fa*
# 	rm ./data/*_repeats.fa*
# fi

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

	# for i in "${SPACER_FILES[@]}"
	# 	do
	# 		echo "Listing spacer files: ${i}"
	# 	done
}

get_repeat_files () {
	REPEAT_FILES=()
	while IFS=  read -r -d $'\0'; do
    	REPEAT_FILES+=("$REPLY")
	done < <(find . -maxdepth 1 -name '*_repeats.fa' -print0)
	echo "Total repeat files available in directory: ${#REPEAT_FILES[@]}"

	# for i in "${REPEAT_FILES[@]}"
	# 	do
	# 		echo "Listing repeat files: ${i}"
	# 	done
}

# split_files () {
# 	fileList=("$@")
# 	for file in "${fileList[@]}"
# 		do
# 			# echo "File from filelist---------->>: ${file}"
# 			item_array=()

# 			awk 'NR % 2 == 1' "${file}" > "${file}.tmp"
# 			# echo "File: ${file}.tmp"
			
# 			while IFS='' read -r line || [[ -n "${line}" ]]; do
#     			# echo "Text read from file: ${line}"
#     			IFS='_ ' read -r -a line_array <<< "${line}"
#     			last_line_item="${line_array[@]: -1:1}"
#     			item_array+=($last_line_item)
# 			done < "${file}.tmp"

# 			# echo "Spacer array: ${item_array[*]}"
# 			newfile_start=1
# 			newfile_count=1
# 			iterator=1
# 			last_idx=0
# 			for idx in "${item_array[@]}"
# 			do
# 				# echo "testing spacer idx | last_idx: ${idx} ${last_idx}"
# 				if [ "$idx" -le "$last_idx" ]; then ## -lt vs -le (<= vs <)
# 					# echo "newstart: $newfile_start"
# 					# echo "iterator: $iterator"
# 					lower_range=$((newfile_start*2-1))
# 					upper_range=$((iterator*2-2))
# 					# echo "creating file from : $lower_range to $upper_range"
# 					awk "NR >= $lower_range && NR <= $upper_range" "${file}" > "${file}.${newfile_count}"
# 					((newfile_count++))
# 					newfile_start=$iterator
# 				fi
# 				last_idx=$idx
# 				((iterator++))
# 			done

# 			lower_range=$((newfile_start*2-1))
# 			upper_range=$((iterator*2-2))
# 			# echo "FINAL: creating file from : $lower_range to $upper_range"
# 			awk "NR >= $lower_range && NR <= $upper_range" "${file}" > "${file}.${newfile_count}"

# 		done
# }

# check_split_files () {
# 	#if split arg, split loci (spacers & repeats) into separate files
# 	if [[ -n $xopt && $xopt -gt 0 ]]; then
# 		echo "checking split files"
# 		get_spacer_files
# 		split_files "${SPACER_FILES[@]}"

# 		get_repeat_files
# 		split_files "${REPEAT_FILES[@]}"

# 		#remove working tmp files
# 		rm *.tmp

# 		if [[ -n $copt && $copt -gt 0 ]]; then
# 			rm *_spacers.fa
# 			rm *_repeats.fa
# 		fi
# 	fi
# }

remove_null_results () {
	if [[ -z ${iopt} ]]; then
		echo "Removing null results..."
		find . -maxdepth 1 -empty -delete
	fi
}


# GENOME_FILES=()
# while IFS=  read -r -d $'\0'; do
#     GENOME_FILES+=("$REPLY")
# done < <(find . -maxdepth 1 -name '*.fasta' -print0)
# echo "Total .fasta files available in directory: ${#GENOME_FILES[@]}"


#check for single file parameter
INPUT_FILE=""
if [ ! -z $fopt ]; then 
	echo "Setting input file: "${fopt}
	INPUT_FILE="$fopt"

	## Copy from R working directory to data directory ##
	#mv "${INPUT_FILE}" "${SOURCE_DIR}"/data/meta.fasta
	#INPUT_FILE="${SOURCE_DIR}"/data/meta.fasta

	${MINCED} -spacers ${MIN_REPEATS}${MIN_REPEAT_LEN}${MAX_REPEAT_LEN}${MIN_SPACER_LEN}${MAX_SPACER_LEN} ${INPUT_FILE} ${INPUT_FILE}.crisprs
	extract_repeats ${INPUT_FILE}.crisprs
	## check_split_files
	remove_null_results
	## execute_spacer_analysis ${INPUT_FILE}_spacers.fa
fi

echo "Minced parsing completed @: " $(date)

exit 1



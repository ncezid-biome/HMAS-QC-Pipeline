#!/bin/bash


# input file format
# 2010K_2370.OG0000318primerGroup3
# 2010K_2370.OG0002922primerGroup0
usage() { echo "Usage: $0 <-c  specify the count_table file (full format required)>" \
				        " <-i specify the input file (with one sample.primer per row)>" 1>&2; exit 1; }

while getopts "c:i:" o; do
	case $o in
		c) COUNT_FILE=${OPTARG} ;;
		i) INPUT_FILE=${OPTARG} ;;
		*) usage ;;
	esac
done

if [[ ! -f "${COUNT_FILE}" || ! -f "${INPUT_FILE}" ]];then
	usage
fi

LINES=$(cat $INPUT_FILE)
for line in $LINES
do
	# check the first line of full format count table file, find the # of field for the matched sample.primer group
	col_num=$(head -n 1 $COUNT_FILE | \
			  awk -v var="$line" '{ for (f=1; f<=NF; f++) { if($f == var) {print f} } }')

	echo "for $line, matched column number in $(basename $COUNT_FILE) is: $col_num"

	if [ -z "$col_num" ]; then
		echo "the total abuandance count is: 0"
	else
		# cut that exact column, remove first row (tail -n +2), count if the value is greater than 0
		cut -f $col_num $COUNT_FILE | tail -n +2 | awk ' $0 > 0 { n++;print "abuandance value is", $0, "at line # ",NR+1  } \
	                                           END { print "the total abuandance count is:", n }'
	fi
done

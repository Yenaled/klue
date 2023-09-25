#!/bin/bash

### Setup ###
klue="./src/klue"
test_dir="./func_tests"

cmdexec() {
	cmd="$1"
	expect_fail="$2"
	printf "$cmd\n"
	if eval "$cmd 2> /dev/null 1> /dev/null"; then
		if [ "$expect_fail" = "1" ]; then
			printf "^[Failed - Returned normally but expected failure]\n"
			exit 1
		else
			printf "^[OK]\n"
		fi
	else
		if [ "$expect_fail" = "1" ]; then
			printf "^[OK - Exited non-zero as expected]\n"
		else
			printf "^[Failed]\n"
			exit 1
		fi
	fi
}

checkcmdoutput() {
	cmd="$1"
	correct_md5="$2"
	printf "$cmd\n"
	output_md5=$(eval "$cmd"|md5sum|awk '{ print $1 }')  # Use 'md5 -r' instead of 'md5sum' on Mac
	if [ "$output_md5" = "$correct_md5" ]; then
		printf "^[Output OK]\n"
	else
		printf "^[Output incorrect! Expected: "
		printf "$correct_md5"
		printf " Actual: "
		printf "$output_md5"
		printf "]\n"
		exit 1
	fi	
}

# Test that program can be run

cmdexec "$klue version"

# Test basic run

checkcmdoutput "$klue distinguish -t 1 -M 1,1 -p $test_dir/test_1.fq.gz $test_dir/test_2.fq.gz|wc -c|tr -d ' '" 96a9f3ee62e50cdc0f6e4afe6fef0ce9
checkcmdoutput "$klue distinguish -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" a8ee4140c073de1ec851bc2b133d5b8b 
checkcmdoutput "$klue distinguish --all -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" a8ee4140c073de1ec851bc2b133d5b8b
checkcmdoutput "$klue distinguish --all-but-one -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" a8ee4140c073de1ec851bc2b133d5b8b
checkcmdoutput "$klue distinguish --combinations -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 897316929176464ebc9ad085f31e7284
# Test extend run

checkcmdoutput "$klue distinguish --extend -t 1 -M 1,1 -p $test_dir/test_1.fq.gz $test_dir/test_2.fq.gz|wc -c|tr -d ' '" e57b5c9f3b511d0e701fe19d51d03f48
checkcmdoutput "$klue distinguish --extend -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 1905d16b771130d54072b7f3db75f98a 
checkcmdoutput "$klue distinguish --extend --all -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 1905d16b771130d54072b7f3db75f98a
checkcmdoutput "$klue distinguish --extend --all-but-one -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 1905d16b771130d54072b7f3db75f98a
checkcmdoutput "$klue distinguish --extend --combinations -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 897316929176464ebc9ad085f31e7284


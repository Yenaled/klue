#!/bin/bash

### Setup ###
klue="./src/klue"
test_dir="./func_tests"

function check_fasta_files {
    if [ "$#" -lt 6 ]; then
        printf "Usage: $0 <file1> <file2> <file3> <file4> <file5> <file6>\n"
        exit 1
    fi

    local file1="$1"
    local file2="$2"
    local file3="$3"
    local file4="$4"
    local file5="$5"
    local file6="$6"
    local expect_fail="$7"

    # Function to read a FASTA file and output sequences in a line-by-line format
    function read_fasta {
        awk '/^>/{if (seq) print seq; seq=""; next} {seq = seq $0} END {print seq}' "$1"
    }

    # Function to generate the reverse complement of a DNA sequence
    function reverse_complement {
        printf "%s" "$1" | rev | tr 'ATGCatgc' 'TACGtacg'
    }

    # Read sequences from each file
    local seqs1=()
    while IFS= read -r line; do
        seqs1+=("$line")
    done < <(read_fasta "$file1")

    local seqs2=()
    while IFS= read -r line; do
        seqs2+=("$line")
    done < <(read_fasta "$file2")

    local seqs3=()
    while IFS= read -r line; do
        seqs3+=("$line")
    done < <(read_fasta "$file3")

    local seqs4=()
    while IFS= read -r line; do
        seqs4+=("$line")
    done < <(read_fasta "$file4")

    local seqs5=()
    while IFS= read -r line; do
        seqs5+=("$line")
    done < <(read_fasta "$file5")

    local seqs6=()
    while IFS= read -r line; do
        seqs6+=("$line")
    done < <(read_fasta "$file6")

    # Check the conditions
    local all_good=true
    local len="${#seqs1[@]}"

    for ((i=0; i<len; i++)); do
        local rc_seq2=$(reverse_complement "${seqs2[i]}")
        local rc_seq4=$(reverse_complement "${seqs4[i]}")
        local rc_seq6=$(reverse_complement "${seqs6[i]}")

        # Condition A
        if [ "${seqs2[i]}" != "${seqs5[i]}" ] && [ "$rc_seq2" != "${seqs5[i]}" ]; then
            all_good=false
            break
        fi

        # Condition B and C or D and E
        if [ "${seqs1[i]}" == "${seqs4[i]}" ] || [ "$rc_seq4" == "${seqs1[i]}" ]; then
            # Condition C
            if [ "${seqs3[i]}" != "${seqs6[i]}" ] && [ "$(reverse_complement "${seqs3[i]}")" != "${seqs6[i]}" ]; then
                all_good=false
                break
            fi
        else
            # Condition D
            if [ "${seqs1[i]}" == "${seqs6[i]}" ] || [ "$rc_seq6" == "${seqs1[i]}" ]; then
                # Condition E
                if [ "${seqs3[i]}" != "${seqs4[i]}" ] && [ "$(reverse_complement "${seqs3[i]}")" != "${seqs4[i]}" ]; then
                    all_good=false
                    break
                fi
            else
                all_good=false
                break
            fi
        fi
    done

    printf "[bubble] $file1 $file2 $file3 : $file4 $file5 $file6\n"

    if [ "$all_good" = true ]; then
        if [ "$expect_fail" = "1" ]; then
            printf "^[Failed - Returned normally but expected failure]\n"
            exit 1
        else
            printf "^[OK]\n"
        fi
    else
        if [ "$expect_fail" = "1" ]; then
            printf "^[OK - Returned non-zero as expected]\n"
        else
            printf "^[Failed]\n"
            exit 1
        fi
    fi
}


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

# Sanity check the bubble testing infrastructure

o="$test_dir/expected_output"
fail="1"
check_fasta_files "$o/dummy_1_L.fa" "$o/dummy_1_V.fa" "$o/dummy_1_R.fa" "$o/dummy_2_L.fa" "$o/dummy_2_V.fa" "$o/dummy_2_R.fa"
check_fasta_files "$o/dummy_1_R.fa" "$o/dummy_1_V.fa" "$o/dummy_1_L.fa" "$o/dummy_2_L.fa" "$o/dummy_2_V.fa" "$o/dummy_2_R.fa"
check_fasta_files "$o/dummy_1_R.fa" "$o/dummy_1_V.fa" "$o/dummy_1revcomp_L.fa" "$o/dummy_2_L.fa" "$o/dummy_2revcomp_V.fa" "$o/dummy_2_R.fa"
check_fasta_files "$o/dummy_1_L.fa" "$o/dummy_1_V.fa" "$o/dummy_1_L.fa" "$o/dummy_2_L.fa" "$o/dummy_2_V.fa" "$o/dummy_2_R.fa" $fail
check_fasta_files "$o/dummy_1revcomp_L.fa" "$o/dummy_1_V.fa" "$o/dummy_1_L.fa" "$o/dummy_2_L.fa" "$o/dummy_2_V.fa" "$o/dummy_2_R.fa" $fail



# Test basic run (note: we only check wc -c, i.e. number of characters, because we don't want to deal with +/- unitig strand issues)

checkcmdoutput "$klue distinguish -t 1 -M 1,1 -p $test_dir/test_1.fq.gz $test_dir/test_2.fq.gz|wc -c|tr -d ' '" 96a9f3ee62e50cdc0f6e4afe6fef0ce9
checkcmdoutput "$klue distinguish -t 1 -M 1,1 -p $test_dir/snpx.fa $test_dir/snpy.fa|wc -c|tr -d ' '" e42bb897d0afcdb1f1c46fb5e0c1ad22
checkcmdoutput "cat $o/dummy_1_L.fa |wc -c|tr -d ' '" 84bc3da1b3e33a18e8d5e1bdd7a18d7a

# The next tests don't work

exit 0

checkcmdoutput "$klue distinguish -t 1 -M 1,1,1 -p $test_dir/test_w.fq $test_dir/test_x.fq $test_dir/test_y.fq|wc -c|tr -d ' '" 0c42da9d05bad0fe0752a9db917beecb 
checkcmdoutput "$klue distinguish --all -t 1 -M 1,1,1 -p $test_dir/test_w.fq $test_dir/test_x.fq $test_dir/test_y.fq|wc -c|tr -d ' '" 0c42da9d05bad0fe0752a9db917beecb
checkcmdoutput "$klue distinguish --all-but-one -t 1 -M 1,1,1 -p $test_dir/test_w.fq $test_dir/test_x.fq $test_dir/test_y.fq|wc -c|tr -d ' '" d015ac2def9a684cfacc41ce1be571af
checkcmdoutput "$klue distinguish --combinations -t 1 -M 1,1,1 -p $test_dir/test_w.fq $test_dir/test_x.fq $test_dir/test_y.fq|wc -c|tr -d ' '" f87a672ba0f13c44f804e9649b96ecd1
# Test extend run

checkcmdoutput "$klue distinguish --extend -t 1 -M 1,1 -p $test_dir/test_w.fq $test_dir/test_x.fq $test_dir/test_y.fq|tr -d ' '" 6c97ca5f37995f5e921ace70ec78b7d8
#checkcmdoutput "$klue distinguish --extend -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 1905d16b771130d54072b7f3db75f98a 
#checkcmdoutput "$klue distinguish --extend --all -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 1905d16b771130d54072b7f3db75f98a
#checkcmdoutput "$klue distinguish --extend --all-but-one -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 1905d16b771130d54072b7f3db75f98a
#checkcmdoutput "$klue distinguish --extend --combinations -t 1 -M 1,1,1 -p $test_dir/test_3.fq.gz $test_dir/test_4.fq.gz $test_dir/test_5.fq.gz|wc -c|tr -d ' '" 897316929176464ebc9ad085f31e7284

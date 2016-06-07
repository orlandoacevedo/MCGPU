#!/bin/bash
#
# timeit.sh -- Times various MCGPU strategies on a couple of sample input files.
#

# Directory containing test files
DIR=resources/exampleFiles
# Test files on which to compute timings
TESTFILES="indole4000.config meoh8788.config"
#TESTFILES="indole267.config meoh500.config"

function run_test {
	name=$1
	args=$2

	printf '%30s' "$name"
	printf '    '
	for file in $TESTFILES; do
		runtime=`bin/metrosim $args "$DIR/$file" | grep 'Run Time:' | sed -e 's/[^0-9.]//g'`
		printf '%20.04f' $runtime
		printf '    '
	done
	printf '\n'
}

echo 'MCGPU will be run on each of the files listed below, once per test.'
echo 'Each run may take several minutes.  Please be patient.'
echo 'Runtimes are in seconds.'
echo ''

# Display headers (config filenames)
printf '%30s' ''
printf '    '
for file in $TESTFILES
do
	printf '%20s' "$file"
	printf '    '
done
printf '\n'

# Display separators (-------)
printf '%30s' ''
printf '    '
for file in $TESTFILES
do
	printf '%20s' '--------------------'
	printf '    '
done
printf '\n'

# Now run tests and display runtimes
run_test 'Parallel, Proximity Matrix:   ' '-p -S proximity-matrix'
run_test 'Parallel, Brute Force:        ' '-p -S brute-force'
run_test 'Serial, Proximity Matrix:     ' '-s -S proximity-matrix'
run_test 'Serial, Brute Force:          ' '-s -S brute-force'

exit 0

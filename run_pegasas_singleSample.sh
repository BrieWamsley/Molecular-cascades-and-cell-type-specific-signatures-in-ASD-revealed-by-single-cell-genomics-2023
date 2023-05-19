#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
for i in `cat sampleList.txt`
	do mkdir $i
	pushd $i
	export $i
	python singlePegasusAnalysis.py
	popd
done

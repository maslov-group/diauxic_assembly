#!/bin/bash
for index in {1..100}; do
	sbatch sbatchjob.pbs "$index"
done

#!/bin/bash
snakemake -npj1 --dag $1 | dot -Tsvg > dag.svg

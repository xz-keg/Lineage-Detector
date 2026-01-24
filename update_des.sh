#!/usr/bin/env bash
set -euo pipefail

curl -L -o des.json "https://nextstrain.org/charon/getDataset?prefix=staging/nextclade/sars-cov-2/"

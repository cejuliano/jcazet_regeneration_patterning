#!/bin/sh

uropa -i resources/full_consensus.json

rm resources/full_consensus_allhits.*
rm resources/full_consensus_finalhits.bed

uropa -i resources/untreated_consensus.json

rm resources/untreated_consensus_allhits.*
rm resources/untreated_consensus_finalhits.bed
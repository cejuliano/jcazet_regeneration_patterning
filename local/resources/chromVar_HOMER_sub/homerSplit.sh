#! /bin/bash

gcsplit -z --digits=3 --quiet --prefix=sub ../chromVar_HOMER.motifs "/>/" "{*}"

for f in sub*; do
    d="$(head -n 1 "$f" | cut -f 2 | cut -d '/' -f 1-2).txt"
    d="${d/\//_}"
    if [ ! -f "$d" ]; then
        mv "$f" "$d"
    fi
done

rm sub*
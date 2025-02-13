#!/bin/bash

# Check if URL is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <URL>"
    exit 1
fi

# Download the webpage
wget -q -O - "$1" | \
# Extract all links
grep -oP '(?<=HREF=")[^"]*' | \
# Filter out only files (not directories)
grep -E '\.[a-zA-Z0-9]{2,4}$' | \
# Download each file
while read -r file; do
    echo file
    wget "http://astro.if.ufrgs.br/evol/evolve/hansen/StellarEvolnDemo/$file"
done



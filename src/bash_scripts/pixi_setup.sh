pixi init

pixi add --platform "osx-arm64" python=3.9
pixi add --platform "linux-64" python=3.9

# Add a channel to the environment
pixi project channel add bioconda

# Maurice to add to Linux machine

pixi add --platform "linux-64" "cutadapt==5.0"

pixi add --platform "linux-64" "seqkit==2.10.0"

pixi add --platform "linux-64" "chopper==0.9.2"

pixi add --platform "linux-64" "grepq"

pixi add --platform "linux-64" "r-optparse"

pixi add --platform "linux-64" bioconductor-decipher

pixi add --platform "linux-64" "picrust2==2.6.2"

# Ellie to add to Mac machine
pixi add --platform "osx-arm64" "bioconductor-decipher"

pixi add --platform "osx-arm64" "picrust2==2.6.2"

# Ellie to remove from Mac machine
pixi remove --platform "osx-arm64" emu

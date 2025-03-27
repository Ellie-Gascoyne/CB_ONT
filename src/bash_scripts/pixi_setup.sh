pixi init

pixi add --platform "osx-arm64" python=3.9
pixi add --platform "linux-64" python=3.9

# Add a channel to the environment
pixi project channel add bioconda

pixi add --platform "linux-64" "cutadapt==5.0"

pixi add --platform "linux-64" "seqkit==2.10.0"

pixi add --platform "linux-64" "chopper==0.9.2"

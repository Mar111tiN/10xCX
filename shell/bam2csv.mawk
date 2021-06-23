#!/bin/sh

# v1.0
# only optional argument is a comma-separated list of additional tags 
# for 10x: CB,CY,UB,UY (default)
# for UMI-tools: RX,QX,MC,MD,RG,NM,AS,XS,MI,cD,cE,cd,ce
 ####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # output depth
        -d|--debug)
        debug=1;
        shift;
        ;;
        # use the following tags
        -t|--tags)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            tags=$2;
            shift 2;
        else
            echo "<bam2csv> Error: tags argument is missing\n[-t|--tags (default=CB,CY,UB,UY)]" >&2
            exit 1;
        fi;;
        # output the following columns
        -c|--cols)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            cols=$2;
            shift 2;
        else
            echo "<bam2csv> Error: cols argument is missing\n[-c|--cols (default=QNAME,FLAG,CHR,START,MAPQ,CIGAR,LENGTH,SEQ,QUAL)]" >&2
            exit 1;
        fi;;
        -*|--*=) # unsupported flags
        echo "<bam2csv> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *) # preserve positional arguments
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done


mawk '
BEGIN {
    ################ BAM2CSV TRANSLATOR ###################

    # define the header columns

    # here go the standard 11 tab-separated fields as described in the SAM-specs
    bamspex="QNAME,FLAG,CHR,START,MAPQ,CIGAR,RNEXT,LENGTH,PNEXT,SEQ,QUAL";
    # get the SPEX array to get the data from
    # SPEX[1] = "QNAME"
    # SPEX[2] = "FLAG"
    # --
    # SPEX[11] = "QUAL"
    spexCount = split(bamspex,SPEX,",");
    # here come the actual output fields. omitting RNEXT and PNEXT
    # can be flexibly changed
    spexCols="'${cols:-QNAME,FLAG,CHR,START,MAPQ,CIGAR,LENGTH,SEQ,QUAL}'";
    spexSelectCount = split(spexCols,XX,",");
    ######### TAGS #####################
    # get extra tags from command argument $1 (with default for 10x)
    tags="'${tags:-CB,CY,UB,UY}'";
    # and paste with the outFields
    cols = spexCols "," tags
    # split the cols string into the DATAFIELDS array
    fieldCount = split(cols,FIELDS,",");
    # FIELDS[1] = "QNAME"
    # ..
    # FIELDS[9] = "QUAL"
    # FIELDS[10] = "CB"
    # ..
    # reverse the array for easy check
    for (i = 0; i++ < fieldCount;) {
      COLS[FIELDS[i]] = i;
    }
    # COLS["QNAME"] = 1 ....
    # COLS["CB"] = 10

    ######## HEADER OUTPUT #############
    for (col = 0; col++ < fieldCount-1;) {printf("%s\t",FIELDS[col])}
    printf("%s\n",FIELDS[fieldCount]);
}

##### PROCESS LINES #################
{
    # get normal fields with COLS check
    for (col = 0; col++ < spexCount;) {
      if (SPEX[col] in COLS) {
        Data[SPEX[col]]=$col;
      }
    }
    # get TAGS 
    for (col=spexSelectCount;col++<fieldCount;) {
      # create the pattern from the DataField and the value format field  (like CB:Z:ATTGCCTG)
      pattern=FIELDS[col] ":[ZiB]:[^\\t]+\\t?[$]?";
      if (match($0, pattern)) {
        Data[FIELDS[col]]=substr($0,RSTART+5,RLENGTH-6);
      } else {
        Data[FIELDS[col]]=".";
      }
    }

    ######## OUTPUT #################
    for (col = 0; col++ < fieldCount - 1;) {
      if (FIELDS[col] in COLS) {
        printf("%s\t",Data[FIELDS[col]]);
      }
    }
    printf("%s\n",Data[FIELDS[fieldCount]]);
}'
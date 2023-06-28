#TOBIAS - calculate footprint scores
#load python/3.6.1
#load samtools
#Yuyan 11-09-19


START=$SECONDS  


INDIR=$1
PEAK=$2
OUTDIR=$3



for bw in $INDIR/*_corrected.bw; do 
    name_prefix=$(basename "$bw")
    name_prefix=${name_prefix/_corrected.bw/}	
    out_bw=$OUTDIR'/'$name_prefix'.footprint.bw'
    
    if [[ ! -f $out_bw ]]; then
    TOBIAS FootprintScores --signal $bw \
    --regions $PEAK \
    --output $out_bw --cores 4
    fi
done

sleep 10
# hack to ensure job lasts 10 minutes to ensure no throttling
END=$SECONDS
ELAPSED=$((END-START))
echo $ELAPSED
if [ $ELAPSED -lt 600 ]; then
  TOSLEEP=$((600 - ELAPSED))
  sleep $TOSLEEP
fi


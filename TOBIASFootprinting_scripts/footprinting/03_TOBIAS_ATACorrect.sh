#TOBIAS - correct Tn5 insertion bias from filtered bam files
#load python/3.6.1
#load samtools
#Yuyan 11-09-19



START=$SECONDS  
INDIR=$1
PEAK=$2
GENOME=$3
OUTDIR=$4


for bam in $INDIR/*.sorted.rmDup.bam; do 
    name_prefix=$(basename "$bam")
    name_prefix=${name_prefix/.bam/}	
    correct_bam=$OUTDIR'/'$name_prefix'_corrected.bw'
    
    if [[ ! -f $correct_bam ]]; then
    
    TOBIAS ATACorrect --bam $bam \
    --genome $GENOME --peaks $PEAK \
    --outdir $OUTDIR --cores 8
    
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


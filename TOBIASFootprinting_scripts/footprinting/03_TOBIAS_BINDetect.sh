#TOBIAS - differential analysis of motif footprints
#load python/3.9.6
#Yuyan 11-09-19


START=$SECONDS  

MOTIF=$1
INDIR=$2
GENOME=$3
PEAK=$4
OUTDIR=$5


TOBIAS BINDetect --motifs $MOTIF \
                 --signals $INDIR'/'*.sorted.rmDup.footprint.bw \
                 --genome $GENOME \
                 --peaks $PEAK \
                 --outdir $OUTDIR --cores 4 --skip-excel

sleep 10
# hack to ensure job lasts 10 minutes to ensure no throttling
END=$SECONDS
ELAPSED=$((END-START))
echo $ELAPSED
if [ $ELAPSED -lt 600 ]; then
  TOSLEEP=$((600 - ELAPSED))
  sleep $TOSLEEP
fi


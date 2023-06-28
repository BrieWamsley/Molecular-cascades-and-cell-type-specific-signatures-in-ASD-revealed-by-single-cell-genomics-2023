#!/bin/bash
## index and collect qc for bam files
#Yuyan Cheng
#11-04-22      

START=$SECONDS 

INDIR=$1
OUTDIR=$2
MAPQ_THR=30


#sort and index
for bam in $INDIR/*.bam;do
	
	name_prefix=$(basename "$bam")
    name_prefix=${name_prefix/.bam/}
    name_prefix=$OUTDIR'/'$name_prefix
    
    # temporary file names
    
    del_mitch_read_bamfile=$name_prefix'.del_Mitch.bam'
    uniq_mapped_read_bamfile=$name_prefix'.UniqMappedRead.bam'
    del_mapq_bam=$name_prefix'.del_MapQ.bam'
    fixmate_bam=$name_prefix'.fixmate.bam'
    sorted_rmdup_bam=$name_prefix'.sorted.rmDup.bam'
    
   
	if [ ! -f $sorted_rmdup_bam ]; then
	  
   	  
	#remove mitonchondria read
	  
	  if [ ! -f $del_mitch_read_bamfile ]; then
	     samtools view -h $bam | sed '/chrM/d;/random/d;/chrUn/d' - | samtools view -Shb -@ 8 -m 4G - > $del_mitch_read_bamfile 	
	  fi
	
	# create BAM file consisting of the uniquely mapped reads
	# the flag 1804 = read unmapped, mate unmapped, not primary alignment, read quality low, optical duplicate  
	  if [ ! -f $uniq_mapped_read_bamfile ]; then
		samtools view -hb -F 1804 $del_mitch_read_bamfile > $uniq_mapped_read_bamfile
		
	  fi
	  
	# remove low-quality reads
	
	  if [[ ! -f $del_mapq_bam ]]; then
	    samtools view -hb -q $MAPQ_THR $uniq_mapped_read_bamfile > $del_mapq_bam
	    
	  fi
	  
	
	  # remove duplicates
   	  if [[ ! -f $sorted_rmdup_bam ]]; then
		# mate score tag missing; need fixmate 
		# requires sort by group name (-n)  	  
   	  	samtools sort -n -@ 8 -m 4G $del_mapq_bam | samtools fixmate -m -@ 8 - $fixmate_bam
   	  	
   	  	# remove duplicate, requires sort by coordiante
	    samtools sort -@ 8 -m 4G $fixmate_bam | samtools markdup -r -@ 8 - $sorted_rmdup_bam
	  
	  
   	  	if [[ ! -f $sorted_rmdup_bam'-m.bai' ]]; then
	    	samtools index $sorted_rmdup_bam
      	elif [[ $sorted_rmdup_bam'.bai' -ot $sorted_rmdup_bam ]]; then
			# here -ot corresponds to "older than"
	    	samtools index $sorted_rmdup_bam
      	fi
    	
	  fi
	  
	  if [[ -f $sorted_rmdup_bam ]]; then
	   
	    rm $fixmate_bam $del_mitch_read_bamfile $uniq_mapped_read_bamfile $del_mapq_bam
	  fi
	  
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




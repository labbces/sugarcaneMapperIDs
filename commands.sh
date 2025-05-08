for i in *_transcript_4120596.fasta.gz;
do
       	echo $i; 
	OUTDIR=${i/.fasta.gz}; 
	if [[ ! -e FCS_outputdir/$OUTDIR ]]; 
	then
	       	mkdir -p FCS_outputdir/$OUTDIR; 
		./run_fcsadaptor.sh --fasta-input $i --output-dir FCS_outputdir/$OUTDIR --euk --container-engine singularity --image fcs-adaptor.sif; 
	fi; 
	if [[ -e FCS_outputdir/$OUTDIR ]];
	then
		if [[ ! -f ${OUTDIR}.fix.fasta.gz ]];
		then
			python3 process_contaminants.py --contaminants FCS_outputdir/$OUTDIR/fcs_adaptor_report.txt --fasta $i --output ${OUTDIR}.fix.fasta.gz --summary ${OUTDIR}.fix.summary.txt --offset 20
		fi
	fi
done

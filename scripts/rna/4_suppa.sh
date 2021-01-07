#!/bin/bash
# 4_suppa.sh
# Ashley Mae Conard
# Last Mod: 7/5/2019
# Purpose: Run various variants of Suppa to to identify differential splicing.

if [ $# -lt 2 ]; then
	echo $0: "Usage: ./4_suppa.sh /PATH/TO/FOLDER/ (e.g. merged_2-4_cf_cm) COL_BEGINNING (e.g. 2-4_cf_cm)"
	exit 1
fi

# Input arguments
FOLDER=$1
COL_BEGINNING=$2

suppa.py psiPerEvent -i /data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/flybase_events/flybase.events.ioe -e $FOLDER/iso_tpm_merged.txt -o $FOLDER/earlyEmbryo_events 2> $FOLDER/errors_iso_tpm_merged.txt

Rscript /data/compbio/aconard/splicing/scripts/SUPPA-2.3/scripts/split_file.R $FOLDER/iso_tpm_merged.txt ${COL_BEGINNING}clampRNAi.1,${COL_BEGINNING}clampRNAi.2,${COL_BEGINNING}clampRNAi.3,${COL_BEGINNING}clampRNAi.4 ${COL_BEGINNING}control.1,${COL_BEGINNING}control.2,${COL_BEGINNING}control.3,${COL_BEGINNING}control.4 clampRNAi_iso.tpm control_iso.tpm -i

Rscript /data/compbio/aconard/splicing/scripts/SUPPA-2.3/scripts/split_file.R $FOLDER/earlyEmbryo_events.psi ${COL_BEGINNING}clampRNAi.1,${COL_BEGINNING}clampRNAi.2,${COL_BEGINNING}clampRNAi.3,${COL_BEGINNING}clampRNAi.4 ${COL_BEGINNING}control.1,${COL_BEGINNING}control.2,${COL_BEGINNING}control.3,${COL_BEGINNING}control.4 clampRNAi_iso_events.psi control_iso_events.psi -i

suppa.py diffSplice -m empirical -gc -i $FOLDER/../flybase_events/flybase.events.ioe -p $FOLDER/clampRNAi_iso_events.psi $FOLDER/control_iso_events.psi -e $FOLDER/clampRNAi_iso.tpm $FOLDER/control_iso.tpm -o $FOLDER/clampRNAi_diffSplice

suppa.py psiPerIsoform -g /data/compbio/aconard/BDGP6/genes.gtf -e $FOLDER/iso_tpm_merged.txt -o $FOLDER/psi_per_fb_trans_isoform 2> $FOLDER/errors_psi_per_fb_trans_isoform.txt

Rscript /data/compbio/aconard/splicing/scripts/SUPPA-2.3/scripts/split_file.R $FOLDER/psi_per_fb_trans_isoform_isoform.psi ${COL_BEGINNING}clampRNAi.1,${COL_BEGINNING}clampRNAi.2,${COL_BEGINNING}clampRNAi.3,${COL_BEGINNING}clampRNAi.4 ${COL_BEGINNING}control.1,${COL_BEGINNING}control.2,${COL_BEGINNING}control.3,${COL_BEGINNING}control.4 clampRNAi_iso_tpm.psi control_iso_tpm.psi -i

suppa.py diffSplice -m empirical -gc -i $FOLDER/../fb_trans.isoforms.ioi -p $FOLDER/clampRNAi_iso_tpm.psi $FOLDER/control_iso_tpm.psi -e $FOLDER/clampRNAi_iso.tpm $FOLDER/control_iso.tpm -o $FOLDER/clampRNAi_diffSplice_iso_dtu
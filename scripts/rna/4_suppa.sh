#!/bin/bash
# 4_suppa.sh
# Ashley Mae Conard
# Last Mod: 7/5/2019
# Purpose: Run two variants of SUPPA to to identify differential splicing (per event and per isoform).

if [ $# -lt 2 ]; then
	echo $0: "Usage: ./4_suppa.sh 
				1) /PATH/TO/RESULTS_FOLDER/ (e.g. 'merged_2-4_cf_cm')
				2) /PATH/TO/FLYBASE_EVENTS_FOLDER/flybase.events.ioe (e.g. 'results/suppa_results_ncbi_trans/flybase_events/flybase.events.ioe')
				3) /PATH/TO/SUPPA/DIR (e.g. '~/Desktop/SUPPA-2.3/')
				4) /PATH/TO/genes.gtf
				5) treatment name (e.g. 'treatment' or 'RNAi')
				6) control name (e.g. 'control' or 'rescue'
				NOTE: NO SPACES OR SPECIAL CHARACTERS IN NAMES besides '-' or '.' ALLOWED."

	exit 1
fi

# Input arguments
RES_FOLDER=$1
FB_EVENTS=$2
SUPPA_DIR=$3
GENES_GTF=$4
TREAT=$5
CONTROL=$6
TRANSISOFORM=/users/pmahable/data/pmahable/data/reference/fb_trans.isoforms.ioi 
# MOD: Added ^^ because of my different file structure. 

# Gather list of treatment and controls
delim=","

# Treatments
# MODIFICATION: Added ${RES_FOLDER} to iso_tpm files. 
head -1 ${RES_FOLDER}/iso_tpm_merged.txt | sed 's/\t/\n/g; /^$/d'| grep $TREAT > header_treatment.txt
head header_treatment.txt
readarray -t treatment_names < header_treatment.txt
echo "Merged ISO TPM experiment replicate names for treatment:${treatment_names[*]}"
rm -f header_treatment.txt
treatment_col_names=$(printf "%s\n$delim\n" "${treatment_names[@]}" | head -n-1 | paste -sd '')

# Controls
head -1 ${RES_FOLDER}/iso_tpm_merged.txt | sed 's/\t/\n/g; /^$/d'| grep $CONTROL > header_control.txt
readarray -t control_names < header_control.txt
echo "Merged ISO TPM experiment replicate names for control:${control_names[*]}"
rm -f header_control.txt
control_col_names=$(printf "%s\n$delim\n" "${control_names[@]}" | head -n-1 | paste -sd '')

# Calculate differential splicing per event
echo ${RES_FOLDER}
echo " suppa.py psiPerEvent --ioe-file $FB_EVENTS --expression-file ${RES_FOLDER}/iso_tpm_merged.txt --output-file ${RES_FOLDER}/splicing_events 2> ${RES_FOLDER}/errors_iso_tpm_merged.txt"
suppa.py psiPerEvent --ioe-file $FB_EVENTS --expression-file ${RES_FOLDER}/iso_tpm_merged.txt --output-file ${RES_FOLDER}/splicing_events 2> ${RES_FOLDER}/errors_iso_tpm_merged.txt
echo "finished 1"
Rscript ${SUPPA_DIR}/scripts/split_file.R ${RES_FOLDER}/iso_tpm_merged.txt $treatment_col_names $control_col_names ${RES_FOLDER}/${TREAT}_iso.tpm ${RES_FOLDER}/${CONTROL}_iso.tpm -i
echo "finished 2"
Rscript ${SUPPA_DIR}/scripts/split_file.R ${RES_FOLDER}/splicing_events.psi $treatment_col_names $control_col_names ${RES_FOLDER}/${TREAT}_iso_events.psi ${RES_FOLDER}/${CONTROL}_iso_events.psi -i
echo "finished 3"
suppa.py diffSplice -m empirical -gc -i $FB_EVENTS -p ${RES_FOLDER}/${TREAT}_iso_events.psi ${RES_FOLDER}/${CONTROL}_iso_events.psi -e ${RES_FOLDER}/${TREAT}_iso.tpm ${RES_FOLDER}/${CONTROL}_iso.tpm -o ${RES_FOLDER}/${TREAT}_vs_${CONTROL}_diffSplice
# MOD: $RES_FOLDER/../flybase_events/flybase.events.ioe replaced with $FB_EVENTS
echo "finished 4" 

# # Calculate differential splicing per isoform
suppa.py psiPerIsoform -g $GENES_GTF -e ${RES_FOLDER}/iso_tpm_merged.txt -o ${RES_FOLDER}/psi_per_fb_trans_isoform 2> ${RES_FOLDER}/errors_psi_per_fb_trans_isoform.txt
echo "finished 5" 
Rscript ${SUPPA_DIR}/scripts/split_file.R ${RES_FOLDER}/psi_per_fb_trans_isoform_isoform.psi $treatment_col_names $control_col_names ${RES_FOLDER}/${TREAT}_iso_tpm.psi ${RES_FOLDER}/${CONTROL}_iso_tpm.psi -i
echo "finished 6" 
echo "suppa.py diffSplice -m empirical -gc -i $TRANSISOFORM -p ${RES_FOLDER}/${TREAT}_iso_tpm.psi ${RES_FOLDER}/$CONTROL_iso_tpm.psi -e ${RES_FOLDER}/${TREAT}_iso.tpm ${RES_FOLDER}/${CONTROL}_iso.tpm -o ${RES_FOLDER}/${TREAT}_vs_${CONTROL}_diffSplice_iso_dtu"
suppa.py diffSplice -m empirical -gc -i $TRANSISOFORM -p ${RES_FOLDER}/${TREAT}_iso_tpm.psi ${RES_FOLDER}/${CONTROL}_iso_tpm.psi -e ${RES_FOLDER}/${TREAT}_iso.tpm ${RES_FOLDER}/${CONTROL}_iso.tpm -o ${RES_FOLDER}/${TREAT}_vs_${CONTROL}_diffSplice_iso_dtu
# MOD: Replacing $RES_FOLDER/../fb_trans.isoforms.ioi with $TRANSISOFORM
echo "finished 7" 

# Full example
# suppa.py psiPerEvent -i /data/compbio/aconard/splicing/results/suppa_results_ncbi_trans/flybase_events/flybase.events.ioe -e $RES_FOLDER/iso_tpm_merged.txt -o $RES_FOLDER/earlyEmbryo_events 2> $RES_FOLDER/errors_iso_tpm_merged.txt
# suppa.py psiPerEvent -i /users/pmahable/data/pmahable/data/flybase_events/flybase.events.ioe -e /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/iso_tpm_merged.txt -o /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/splicing_events 2> /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/errors_iso_tpm_merged.txt

# Rscript /data/compbio/aconard/splicing/scripts/SUPPA-2.3/scripts/split_file.R $RES_FOLDER/iso_tpm_merged.txt ${COL_BEGINNING}clampRNAi.1,${COL_BEGINNING}clampRNAi.2,${COL_BEGINNING}clampRNAi.3,${COL_BEGINNING}clampRNAi.4 ${COL_BEGINNING}control.1,${COL_BEGINNING}control.2,${COL_BEGINNING}control.3,${COL_BEGINNING}control.4 clampRNAi_iso.tpm control_iso.tpm -i
# Rscript /users/pmahable/data/pmahable/SUPPA/scripts/split_file.R /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/iso_tpm_merged.txt bf_S2_GFP.1,bf_S2_GFP.2,bf_S2_GFP.3 bm_KC_GFP.1,bm_KC_GFP.2,bm_KC_GFP.3 bf_S2_GFP_iso.tpm bm_KC_GFP_iso.tpm -i

# Rscript /data/compbio/aconard/splicing/scripts/SUPPA-2.3/scripts/split_file.R $RES_FOLDER/earlyEmbryo_events.psi ${COL_BEGINNING}clampRNAi.1,${COL_BEGINNING}clampRNAi.2,${COL_BEGINNING}clampRNAi.3,${COL_BEGINNING}clampRNAi.4 ${COL_BEGINNING}control.1,${COL_BEGINNING}control.2,${COL_BEGINNING}control.3,${COL_BEGINNING}control.4 clampRNAi_iso_events.psi control_iso_events.psi -i
# Rscript /users/pmahable/data/pmahable/SUPPA/scripts/split_file.R /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/splicing_events.psi bf_S2_GFP.1,bf_S2_GFP.2,bf_S2_GFP.3 bm_KC_GFP.1,bm_KC_GFP.2,bm_KC_GFP.3 bf_S2_GFP_iso_events.psi bm_KC_GFP_iso_events.psi -i

# suppa.py diffSplice -m empirical -gc -i $RES_FOLDER/../flybase_events/flybase.events.ioe -p $RES_FOLDER/clampRNAi_iso_events.psi $RES_FOLDER/control_iso_events.psi -e $RES_FOLDER/clampRNAi_iso.tpm $RES_FOLDER/control_iso.tpm -o $RES_FOLDER/clampRNAi_diffSplice
# python3 /users/pmahable/data/pmahable/SUPPA/suppa.py diffSplice -m empirical -gc -i /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/../flybase_events/flybase.events.ioe -p /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/bf_S2_GFP_iso_events.psi /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/.psi -e /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/.tpm /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/bm_KC_GFP_iso.tpm -o /users/pmahable/data/pmahable/cell_line_t2s/time2splice/results/analysis/suppa/merged_gfp_fm/bf_S2_GFP_vs_

# suppa.py psiPerIsoform -g /data/compbio/aconard/BDGP6/genes.gtf -e $RES_FOLDER/iso_tpm_merged.txt -o $RES_FOLDER/psi_per_fb_trans_isoform 2> $RES_FOLDER/errors_psi_per_fb_trans_isoform.txt

# Rscript /data/compbio/aconard/splicing/scripts/SUPPA-2.3/scripts/split_file.R $RES_FOLDER/psi_per_fb_trans_isoform_isoform.psi ${COL_BEGINNING}clampRNAi.1,${COL_BEGINNING}clampRNAi.2,${COL_BEGINNING}clampRNAi.3,${COL_BEGINNING}clampRNAi.4 ${COL_BEGINNING}control.1,${COL_BEGINNING}control.2,${COL_BEGINNING}control.3,${COL_BEGINNING}control.4 clampRNAi_iso_tpm.psi control_iso_tpm.psi -i

# suppa.py diffSplice -m empirical -gc -i $RES_FOLDER/../fb_trans.isoforms.ioi -p $RES_FOLDER/clampRNAi_iso_tpm.psi $RES_FOLDER/control_iso_tpm.psi -e $RES_FOLDER/clampRNAi_iso.tpm $RES_FOLDER/control_iso.tpm -o $RES_FOLDER/clampRNAi_diffSplice_iso_dtu
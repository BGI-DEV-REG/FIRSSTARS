#!/bin/bash
# **********************************************************
# * Author        : HuangFubaoqian
# * Email         : huangbaoqian@genomics.cn
# * Create time   : 2022-08-31 10:31
# * Filename      : pyscenic_from_loom.sh
# * Description   :
# **********************************************************

source /hwfssz5/ST_SUPERCELLS/P21Z10200N0171/USER/zengxiaoqi/software/miniconda3/bin/activate /hwfssz5/ST_SUPERCELLS/P21Z10200N0171/USER/wangyijin/00.Software/miniconda/envs/scenicplus
#source activate /hwfssz5/ST_SUPERCELLS/P21Z10200N0171/USER/fengweimin/software/conda/miniconda3/envs/scenicplus

#default value
input_loom=out.loom
n_workers=20
output_path="output_directory/"

#help function
function usage() {
echo -e "OPTIONS:\n-i|--input_loom:\t input loom file"
echo -e "-n|--n_workers:\t working core number"
echo -e "-o|--output_path:\t output directory path"
echo -e "-h|--help:\t Usage information"
exit 1
}

#get value
while getopts :i:n:o:h opt
do
    case "$opt" in
        i) input_loom="$OPTARG" ;;
        n) n_workers="$OPTARG" ;;
        o) output_path="$OPTARG" ;;
        h) usage ;;
        :) echo "This option -$OPTARG requires an argument."
           exit 1 ;;
        ?) echo "-$OPTARG is not an option"
           exit 2 ;;
    esac
done

database='/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/USER/wangyijin/00.package/03.SCENIC/00.databases/01.mouse'
tfs=${database}/allTFs_mm.txt
feather=${database}/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=${database}/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
pyscenic=pyscenic

# grn
$pyscenic grn \
--num_workers $n_workers \
--output ${output_path}/grn.tsv \
--method grnboost2 \
$input_loom  $tfs

# cistarget
$pyscenic ctx \
${output_path}/grn.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output ${output_path}/ctx.csv \
--num_workers $n_workers   \
--mask_dropouts

# AUCell
$pyscenic aucell \
$input_loom \
${output_path}/ctx.csv \
--output ${output_path}/aucell.loom \
--num_workers $n_workers


#!/bin/bash

# Default values
genomesize=""
threads=20
fastq=""
steps=()
motif="AACCCT"
assembler="peregrine" # Default assembler

# Function to display usage
usage() {
  echo "Usage: $0 -g <genomesize> -t <threads> --fastq <fastq> -m <motif> [-s <steps>] [--choose]"
  echo ""
  echo "Options:"
  echo "  -g    Genome size (required) <number>[g|m|k]"
  echo "  -t    Number of threads (required)"
  echo "  --fastq    Path to the FASTQ file (required)"
  echo "  -m    Telomere motif (required)"
  echo "  -s    Steps to run (optional, default: all steps)."
  echo "        You can specify individual steps or ranges (e.g., 1,2,3 or 1-6 or 2-4)."
  echo "  --choose  Prompt to choose an assembler for the final merge (optional, default: peregrine)."
  echo ""
  echo "Steps:"
  echo "  1. Assembly of the genome using HiCanu"
  echo "  2. Assembly of the genome using NextDenovo"
  echo "  3. Assembly of the genome using Peregrine"
  echo "  4. Assembly of the genome using IPA"
  echo "  5. Assembly of the genome using Flye"
  echo "  6. Assembly of the genome using Hifiasm"
  echo "  7. Copy all assemblies"
  echo "  8. Extract telomere-containing contigs"
  echo "  9. Merge all assemblies"
  echo "  10. Run quast.py for all assembler results"
  echo "  11. Final merge using selected assembler"
  echo "  12. BUSCO analysis"
  echo "  13. Telomere analysis"
  echo ""
  echo "Example:"
  echo "  $0 -g 2g -t 16 --fastq /path/to/reads.fastq -m AACCCT -s 1,3-5 --choose"
  exit 1
}

# Function to expand ranges in the step input (e.g., "1,3-5" becomes "1 3 4 5")
expand_steps() {
  local IFS=','; read -ra ranges <<< "$1"
  for range in "${ranges[@]}"; do
    if [[ $range == *"-"* ]]; then
      local start=${range%-*}
      local end=${range#*-}
      for ((i=start; i<=end; i++)); do
        steps+=("$i")
      done
    else
      steps+=("$range")
    fi
  done
}

# Parse command-line options
while getopts ":g:t:s:m:-:" opt; do
  case $opt in
    g) genomesize="$OPTARG"
    ;;
    t) threads="$OPTARG"
    ;;
    s) expand_steps "$OPTARG"
    ;;
    m) motif="$OPTARG"
    ;;
    -) case "${OPTARG}" in
         fastq) fastq="${!OPTIND}"; OPTIND=$((OPTIND + 1))
         ;;
         choose) 
           echo "Please enter the assembler you want to use for the final merge (caun, nextDenovo, peregrine, ipa, flye, RAFT-hifiasm):"
           read assembler
           ;;
         *) echo "Unknown option --${OPTARG}" >&2; usage
         ;;
       esac
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        usage
    ;;
    :) echo "Option -$OPTARG requires an argument." >&2
       usage
    ;;
  esac
done

# Ensure that required options are set
if [[ -z "$genomesize" || -z "$threads" || -z "$fastq" || -z "$motif" ]]; then
  usage
fi

project=$(basename "$fastq" .fastq)

# Modify run.cfg with genomesize and threads
cat <<EOT >> ./run_${project}.cfg
[General]
job_type = local # local, slurm, sge, pbs, lsf
job_prefix = nextDenovo
task = all # all, correct, assemble
rewrite = yes # yes/no
deltmp = yes
parallel_jobs = ${threads} # number of tasks used to run in parallel
input_type = raw # raw, corrected
read_type = hifi # clr, ont, hifi
input_fofn = ./input_${project}.fofn
workdir = NextDenovo

[correct_option]
read_cutoff = 1k
genome_size = ${genomesize} # estimated genome size
sort_options = -m 20g -t ${threads}
minimap2_options_raw = -t ${threads}
pa_correction = 3 # number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage.
correction_options = -p ${threads}

[assemble_option]
minimap2_options_cns = -t ${threads}
nextgraph_options = -a 1

# see https://nextdenovo.readthedocs.io/en/latest/OPTION.html for a detailed introduction about all the parameters"
EOT

echo $PWD/${fastq} > input_${project}.fofn
echo $PWD/${fastq} > reads_${project}.lst

# Function to check the success of the previous command
check_command() {
  if [ $? -ne 0 ]; then
    echo "Error: Command failed. Exiting."
    exit 1
  fi
}

echo "Activating assembly environment..."
eval "$(conda shell.bash hook)"
conda activate pacbiohifi
check_command

# If no specific steps provided, default to running all steps (1-13)
if [ ${#steps[@]} -eq 0 ]; then
  steps=(1 2 3 4 5 6 7 8 9 10 11 12 13)
fi

# Execute the specified steps
for step in "${steps[@]}"; do
  case $step in
    1)
      echo "Step 1 - Assembly of the genome using HiCanu"
      canu -p canu -d hicanu genomeSize=$genomesize maxThreads=$threads -pacbio-hifi $fastq
      check_command
      ;;
    2)
      echo "Step 2 - Assembly of the genome using NextDenovo"
      nextDenovo run_${project}.cfg
      check_command
      ;;
    3)
      echo "Step 3 - Assembly of the genome using Peregrine"
      pg_asm reads_${project}.lst peregrine-2021
      check_command
      ;;
    4)
      echo "Step 4 - Assembly of the genome using IPA"
      ipa local --nthreads $threads --njobs 1 --run-dir ipa -i $fastq
      check_command
      ;;
    5)
      echo "Step 5 - Assembly of the genome using Flye"
      flye --pacbio-hifi $fastq --out-dir flye --threads $threads
      check_command
      ;;
    6)
      echo "Step 6 - Assembly of the genome using Hifiasm"
      mkdir hifiasm
      cd hifiasm
      hifiasm -o errorcorrect -t$threads --write-ec ../$fastq 2> errorcorrect.log
      check_command
      COVERAGE=$(grep "homozygous" errorcorrect.log | tail -1 | awk '{print $6}')
      hifiasm -o getOverlaps -t$threads --dbg-ovec errorcorrect.ec.fa 2> getOverlaps.log
      check_command
      cat getOverlaps.0.ovlp.paf getOverlaps.1.ovlp.paf > overlaps.paf
      ~/opt/RAFT/raft -e ${COVERAGE} -o fragmented errorcorrect.ec.fa overlaps.paf
      hifiasm -o finalasm -t$threads -r1 fragmented.reads.fasta 2> finalasm.log
      check_command
      awk '/^S/{print ">"$2;print $3}' finalasm.bp.hap1.p_ctg.gfa > RAFT-hifiasm.fasta
      cd ..
      ;;
    7)
      echo "Step 7 - Copy all assemblies"
      cp ./hicanu/canu.contigs.fasta ./canu.result.fasta
      cp ./NextDenovo/03.ctg_graph/nd.asm.fasta ./nextDenovo.result.fasta
      cp ./peregrine-2021/asm_ctgs_m_p.fa ./peregrine.result.fasta
      cp ./ipa/assembly-results/final.p_ctg.fasta ./ipa.result.fasta
      cp ./flye/assembly.fasta ./flye.result.fasta
      cp ./hifiasm/RAFT-hifiasm.fasta ./RAFT-hifiasm.result.fasta
      ;;
    8)
      echo "Step 8 - Extract telomere-containing contigs"
      for fasta in *.result.fasta; do
        seqtk telo -s 1 -m "$motif" $fasta > "${fasta%.result.fasta}.telo.list"
        bash ./extract_contig_T_V3.sh -i "${fasta%.result.fasta}.telo.list" -f $fasta -o "${fasta%.result.fasta}.telo.fasta"
        check_command
      done
      ;;
    9)
      echo "Step 9 - Merge all telo assemblies"
      fasta_files=(*.telo.fasta)
    for file1 in "${fasta_files[@]}"; do
        for file2 in "${fasta_files[@]}"; do
          base1=$(basename "$file1" .telo.fasta)
          base2=$(basename "$file2" .telo.fasta)
          merge_wrapper.py -l 1000000 "$file1" "$file2" --prefix merged_"$base1"_"$base2"
          check_command
        done
      done
      cat merged_*.fasta > allmerged_telo.fasta
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate funannotate
      funannotate sort -i allmerged_telo.fasta -b contig -o allmerged_telo_sort.fasta --minlen 500
      seqtk telo -s 1 -m "$motif" allmerged_telo_sort.fasta > allmerged.telo.list
      bash ./t2t_list.sh -i allmerged.telo.list -o t2t.list
      ~/opt/scripts/faSomeRecords allmerged_telo_sort.fasta t2t.list t2t.fasta
      funannotate clean -i  t2t.fasta -p 30 -o  t2t_clean.fasta --exhaustive
      check_command
      ;;
    10)
      echo "Step 10 - Run quast.py for all assembler results"
      echo "Activating assembly environment..."
      eval "$(conda shell.bash hook)"
      conda activate pacbiohifi
      quast.py *.result.fasta --threads $threads
      check_command
      ;;
    11)
      # Step 11 - Final merge using selected assembler
      echo "Step 11 - Final merge"
      # Add assembler choice logic
      if [[ "$1" == "--choose" ]]; then
          echo "Please enter the assembler you want to use for the final merge (caun, nextDenovo, peregrine, ipa, flye, RAFT-hifiasm):"
          read assembler
      if [[ ! "$assembler" =~ ^(caun|nextDenovo|peregrine|ipa|flye|RAFT-hifiasm)$ ]]; then
        echo "Invalid assembler selected. Exiting."
        exit 1
      fi
      merge_wrapper.py -l 1000000 ${assembler}.result.fasta t2t_clean.fasta --prefix "$assembler"
      check_command
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate funannotate
      funannotate clean -i merged_${assembler}.fasta -p 30 -o merged_${assembler}_clean.fa --exhaustive
      check_command
      funannotate sort -i merged_${assembler}_clean.fa -b contig -o merged_${assembler}_sort.fa --minlen 500
      check_command
      ;;
    12)
      # BUSCO analysis
      echo "Step 12 - BUSCO analysis"
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate busco
      busco -o merged -i merged_${assembler}_sort.fa -l ascomycota_odb10 --cpu $threads -m genome
      check_command
      ;;
    13)
      # Telomere analysis
       echo "Step 13 - Telomere analysis"
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate quartet
      ~/opt/telomere_analysis.sh merged_${assembler}_sort.fa
      ;;
    *)
      echo "Invalid step: $step"
      exit 1
      ;;
  esac
done

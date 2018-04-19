# Allow over-ride
if [ -z "${CONTAINER_IMAGE}" ]
then
    version=$(cat ./_util/VERSION)
    CONTAINER_IMAGE="sd2e/psap:$version"
fi
. _util/container_exec.sh

log(){
    mesg "INFO" $@
}

die() {
    mesg "ERROR" $@
    # AGAVE_JOB_CALLBACK_FAILURE is macro populated at runtime
    # by the platform and gives us an eject button from
    # anywhere the application is running
    ${AGAVE_JOB_CALLBACK_FAILURE}
}

mesg() {
    lvl=$1
    shift
    message=$@
    echo "[$lvl] $(utc_date) - $message"
}

utc_date() {
    echo $(date -u +"%Y-%m-%dT%H:%M:%SZ")
}

# Double check existence of fastq_dir before undertaking
# expensive processing and/or analysis. Fail if not found. 
if [ ! -d "${fastq_dir}" ];
then
    die "inputData ${fastq_dir} not found or was inaccessible"
fi

# Add contents of some child directories to .agave.archive
# Why? Because we don't need to copy the assay and controls
# back out at the end. 
# Agave uses .agave.archive to mark files that were present 
# before the core application logic began running. It's 
# generated automatically when the dependencies are staged
# into place on the executionHost and we're just extending
# that process a bit
inputDir="${fastq_dir}"
for FN in inputDir
do
    echo "${FN}" >> .agave.archive
done

# We have not failed yet. Systems are probably nominal.
# Kick off the analysis
container_exec ${CONTAINER_IMAGE} python /prot_stab_analysis_pipeline/scripts/compute_ec50_values_from_deep_sequencing_data.py --designed_sequences_file ${designed_sequences_file} --fastq_dir ${fastq_dir}  --experimental_summary_file ${experimental_summary_file} --output_dir ${output_dir} --pear_path ${pear_path}

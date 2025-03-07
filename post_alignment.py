import argparse
import subprocess
import time

## Chang input directories 
chrM_file = "path/to/chrM.bed"
picard = "path/to/picard.jar"

def alignment_adjustment(input_bam:str, sample_name: str, output_dir:str) -> None:

    # remove mitochondrial reads
    ## check the percentage of mt reads

    mt_reads_command = f"samtools idxstats {input_bam} | grep 'chrM' | cut -f 3"
    mt_reads_output = subprocess.check_output(mt_reads_command, shell=True).strip()
    mt_reads = int(mt_reads_output) if mt_reads_output else 0

    total_reads_command = f"samtools idxstats {input_bam} | awk '{{SUM += $3}} END {{print SUM}}'"
    total_reads_output = subprocess.check_output(total_reads_command, shell=True).strip()
    total_reads = int(total_reads_output) if total_reads_output else 0

    mt_dna_content = (mt_reads / total_reads) * 100
    content = f'Mitochondrial Reads: {mt_reads}\nTotal Reads: {total_reads}\nPercentage of mtDNA: {mt_dna_content}%\n'
    ## output the percentage statsitics
    subprocess.run(
         f'echo "{content}" > {output_dir}/{sample_name}_chrM_statistics.txt',
         shell=True, check=True
    )
    ## create new bam file with chrM filtered
    cmd_mt_removal = f'bedtools intersect -v -abam {input_bam} -b {chrM_file} > {output_dir}/{sample_name}_noMT.bam'
    subprocess.run(cmd_mt_removal, shell=True, check=True)
    new_bam = f'{output_dir}/{sample_name}_noMT.bam'
    ## reindex the new bam file
    subprocess.run(f"samtools index {new_bam}", shell=True, check=True)

    # add RG group
    cmd_add_rg = f'java -jar {picard} AddOrReplaceReadGroups -I {new_bam} -O {output_dir}/{sample_name}_noMT_rg.bam -LB lib1 -PL Illumina -SM atac -PU unit1'
    subprocess.run(cmd_add_rg, shell=True, check=True)
    new_bam = f'{output_dir}/{sample_name}_noMT_rg.bam'

    # remove duplicate reads
    cmd_dups_removal = f'java -jar {picard} MarkDuplicates -I {new_bam} -O {output_dir}/{sample_name}_noMT_noDups.bam -M {sample_name}_dups.txt -REMOVE_DUPLICATES true'
    subprocess.run(cmd_dups_removal, shell=True, check=True)

    # remove non-unique alignments (set minimum mapping score as 10)
    cmd_nonunique_removal = f'samtools view -b  -q 10  {new_bam}  >  {output_dir}/{sample_name}_noMT_noDups_unique.bam'
    subprocess.run(cmd_nonunique_removal, shell=True, check=True)

    subprocess.run('echo Post-alignment DONE', shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description='Post Alignment Adjustment')
    parser.add_argument('-b', '--bam', type=str)
    parser.add_argument('-s', '--sample_name', type=str)
    parser.add_argument('-o', '--output', type=str)
    args = parser.parse_args()

    alignment_adjustment(args.bam, args.sample_name, args.output)

if __name__ == "__main__":
       start_time = time.time()
       main()
       print('------- %s seconds ---------' % (time.time() - start_time))

import os
import sys


# result_fold = '/slurm/home/zju/zhanglab/linfuqiang/project/03Denmark.birds/02.busco/results/'
result_fold = sys.argv[1]
result_dic = {}

for fold in os.listdir(result_fold):
    if not fold.endswith('.log') and not fold.endswith('.out'):
        sample = fold
        for file in os.listdir(result_fold + fold):
            if file.endswith('.txt'):
                with open(result_fold + fold + '/short_summary.specific.aves_odb10.' + sample + '.txt', 'r') as fin:
                    for line in fin.readlines():
                        if '	C:' in line:
                            single = line.split(':')[2].split('%')[0]
                            duplication = line.split(':')[3].split('%')[0]
                            fragmented = line.split(':')[4].split('%')[0]
                            missing = line.split(':')[5].split('%')[0]
                            result_dic[sample] = [single, duplication, fragmented, missing]
                        if 'Total length' in line:
                            genome_size = line.split('\t')[1]
                            result_dic[sample].append(genome_size)
                        if 'Contigs N50' in line:
                            if line.split()[1] == 'KB':
                                contigs_N50 = str(1024 * int(line.split()[0]))
                            elif line.split()[1] == 'MB':
                                contigs_N50 = str(1024 * 1024 * int(line.split()[0]))
                            result_dic[sample].append(contigs_N50)
                        if 'Scaffold N50' in line:
                            if line.split()[1] == 'KB':
                                scaffold_N50 = str(1024 * int(line.split()[0]))
                            elif line.split()[1] == 'MB':
                                scaffold_N50 = str(1024 * 1024 * int(line.split()[0]))
                            result_dic[sample].append(scaffold_N50)

with open("./busco.result.csv", 'w') as fout:
    fout.write(f"sample,Contig N50(bp),Scaffold N50(bp),Genome size(bp),Single(%),Duplication(%),Fragmented(%),"
               f"Missing(%)\n")
    for key, value in result_dic.items():
        sample_name = key
        fout.write(f"{sample_name},{value[6]},{value[5]},{value[4]},{value[0]},{value[1]},{value[2]},{value[3]}\n")

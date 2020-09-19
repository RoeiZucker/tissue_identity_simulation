import numpy as np
import pandas as pd
from random import shuffle
import ast

def load_data(atlas_path, cpg_snp_path):
    atlas = pd.read_csv(atlas_path)
    names = atlas['0']
    atlas.drop(['0'], axis=1, inplace=True)
    atlas = atlas.rename(index=names)

    sites = pd.read_csv(cpg_snp_path)
    sites.columns = ['cpg_name', 'cpg_address', 'snp_name', 'snp_address', 'snp_freq', 'major', 'minor']
    names = sites['cpg_name']
    sites = sites.rename(index=names)
    sites.drop(['cpg_name'], axis=1, inplace=True)

    # print(atlas)
    # print(sites)

    return atlas, sites



def create_samples(depth, atlas, sites, *args, **kargs):
    """
    :param depth: Sequencing depth
    :param args: Dictionary, where key is tissue name and key is percentage of tissue in mixed sample
    :return:
    """

    separated_samples = {}
    persons_profiles = {}
    for tissue in kargs.keys():
        one_person_sample = {}
        profile = {}
        for cpg in sites.index:
            major = sites.loc[cpg]['major']
            minor = sites.loc[cpg]['minor']
            freq = sites.loc[cpg]['snp_freq']
            snp = np.random.choice([major, minor], p=[freq, 1 - freq])
            one_person_sample[cpg] = []
            profile[cpg] = snp
            C_num = int(depth * atlas.loc[cpg][tissue])
            for i in range(C_num):
                one_person_sample[cpg].append(['C', snp])
            for i in range(depth - C_num):
                one_person_sample[cpg].append(['T', snp])
        separated_samples[tissue] = one_person_sample
        persons_profiles[tissue] = profile
    return separated_samples, persons_profiles


def test(atlas_path, cpg_snp_path, tissue_origin):

    with open(tissue_origin, mode='r') as cpgFile:
        for line in cpgFile:
            tissue_origin_dict  = ast.literal_eval(line)
    atlas, sites = load_data(atlas_path, cpg_snp_path)
    sep_samples, persons_profiles = create_samples(30, atlas, sites, **{'Saliva_edited': 0.33, 'Whole_blood': 0.33, 'Skin':0.33})
    mixed_samples = mix(sep_samples)
    mixed_samples = drop_redundant(mixed_samples)
    detected_tissues = detect(atlas, mixed_samples, persons_profiles)
    # matched_results = match(atlas, detected_tissues, persons_profiles)
    keys = []
    for key in detected_tissues.keys():
        keys.append(key)
    miss_cout = 0
    for a in keys:
        if a not in tissue_origin_dict:
            miss_cout+=1
            continue
        if detected_tissues[a][0] !=tissue_origin_dict[a]:
            del detected_tissues[a]
    # print("dropped count: " + str(miss_cout))

    match(atlas, detected_tissues, persons_profiles)


def mix(samples):
    mixed_samples = {}
    for tissue in samples.keys():
        for site in samples[tissue].keys():
            mixed_samples[site] = mixed_samples.get(site, []) + samples[tissue][site]
    for site in mixed_samples.keys():
        shuffle(mixed_samples[site])
    return mixed_samples


def drop_redundant(samples):
    # print(samples)
    # print("Before: ", len(samples))
    keys_to_drop = []
    for site in samples:
        snp = samples[site][0][1]
        if all(read[1] == snp for read in samples[site]):
            keys_to_drop.append(site)
    for key in keys_to_drop:
        samples.pop(key)
    return samples
    # print("To drop: ", len(keys_to_drop))
    # print("After: ", len(samples.keys()))


def detect(atlas, samples, persons_profiles):
    detected = {}
    for cpg in samples.keys():
        indicative_reads = []
        reads_dict = {}
        for read in samples[cpg]:
            reads_dict[read[1]] = reads_dict.get(read[1], []) + [read]

        for value in reads_dict.values():
            if len(value) < len(indicative_reads) or len(indicative_reads) == 0:
                indicative_reads = value

        methylated = 0.0
        for read in indicative_reads:
            if read[0] == 'C':
                methylated += 1
        methylated_ratio = methylated / len(indicative_reads)
        diff = 1
        detected_tissue = None
        for tissue in atlas.columns:
            if abs(methylated_ratio - atlas.loc[cpg][tissue]) < diff:
                detected_tissue = tissue
                diff = abs(methylated_ratio - atlas.loc[cpg][tissue])
        detected[cpg] = [detected_tissue, indicative_reads[0][1]]
    return detected


def match(atlas, detected_tissues, persons_profiles):
    # result = {}
    global all_val, not_org_val, miss_val, correct_find_dict, incorrect_find_dict


    counter1 = 0
    counter2 = 0
    all = 0
    for cpg in detected_tissues.keys():
        all +=1
        tissue = detected_tissues[cpg][0]
        if persons_profiles.get(tissue, None) == None:
            # print("Tissue that wasn't originally in samples")
            incorrect_find_dict[tissue] += 1
            counter1 += 1
        elif detected_tissues[cpg][1] != persons_profiles[tissue][cpg]:
            counter2 += 1
            incorrect_find_dict[tissue]+=1
            # print("Missdetected")
        else:
            correct_find_dict[tissue]+=1
    # print("All: ", all)
    all_val+=all
    # print("Not org", counter1)
    not_org_val+=counter1
    # print("Miss:", counter2)
    miss_val+=counter2
    # print("correctly diagnosed tissues: ", end='')
    # print(correct_find_dict)
    # print("incorrectly diagnosed tissues: ", end='')
    # print(incorrect_find_dict)

# load_data("/Users/a1/Documents/My Documents/3rd Year Project/Simulation/atlas_06_05_20.tsv",
#           "/Users/a1/Documents/My Documents/3rd Year Project/Simulation/cpg_snp_500_full.csv")

miss_val = 0
not_org_val = 0
all_val = 0
correct_find_dict = {'Urine': 0, 'Sperm': 0, 'Skin': 0, 'Saliva_edited': 0, 'Whole_blood': 0}
incorrect_find_dict = {'Urine': 0, 'Sperm': 0, 'Skin': 0, 'Saliva_edited': 0, 'Whole_blood': 0}
path = "C:\\Users\\roeis\\Dropbox\\study\\lab\\data_files\\identity_tissue_detection"
for i in range(100):
    test("C:\\Users\\roeis\\Dropbox\\study\\lab\\data_files\\raw_data\\current_avrg_matrix.csv", path + "\\cpg_snp_500_full_16_09_2020.csv", path+"\\tissues_orig_16_09_2020.txt")

print("miss val " + str(miss_val))
print("not_org  val " + str(not_org_val))
print("all val " + str(all_val))
print("correctly diagnosed tissues: ", end='')
print(correct_find_dict)
print("incorrectly diagnosed tissues: ", end='')
print(incorrect_find_dict)

import os
import argparse
from mimic3benchmark.subject import read_stays, read_diagnoses
from mimic3benchmark.util import *
import yaml
import collections
import random

table_name_map= {"drgcodes" : "drg", "prescriptions" : "pres", "procedures": "proc", "labevents":"labt", "microbiologyevents":"micro"}
file_header = ['NODE1', 'NODE2', 'COUNT', 'HOURS', 'TYPE']


def write_bio_edges(age, gender,race, pati_abrv, output_dir):
    age = "O" if age > 50 else "Y"
    bio_df = pd.DataFrame(columns=file_header)
    bio_df.loc[len(bio_df)] = [pati_abrv, "age." + age, "1", "", "pati.bio"]
    bio_df.loc[len(bio_df)] = [pati_abrv, "gender." + str(gender), 1, "", "pati.bio"]
    bio_df.loc[len(bio_df)] = [pati_abrv, "eth." + str(race), 1, "", "pati.bio"]

    bio_df.to_csv(os.path.join(output_dir, "Edge.csv"), index_label='NODE1', mode="a", header=False, index=False)

def create_pheno_mappings(args):
    with open(args.phenotype_definitions) as definitions_file:
        definitions = yaml.load(definitions_file)

    code_to_group = {}
    for group in definitions:
        codes = definitions[group]['codes']
        for code in codes:
            if code not in code_to_group:
                code_to_group[code] = group
            else:
                assert code_to_group[code] == group

    id_to_group = sorted(definitions.keys())
    group_to_id = dict((x, i) for (i, x) in enumerate(id_to_group))
    return id_to_group, group_to_id


def write_seq_event(ts_lines, output_dir, pati_abrv, table_abrv_name):

    normal_res = []
    abnormal_res = []
    header = ts_lines[0]
    ts_lines = ts_lines[1:]
    event_times = header.split(',')[1:]

    for line in ts_lines:
        normal_labs = []
        abnormal_labs = []
        count_normal = 0
        count_abnormal = 0
        line_parts = line.split(",")
        lab_id = line_parts[0]
        for lab_ind, lab_val in enumerate(line_parts[1:]):
            if lab_val == "ab":
                abnormal_labs.append(float(event_times[lab_ind]))
                count_abnormal += 1
            elif lab_val == "n":
                try:
                    normal_labs.append(float(event_times[lab_ind]))
                    count_normal += 1
                except:
                    print ("exception in reading a lab value")
                    float(event_times[lab_ind].replace('\U00002013', '-'))

        if len(normal_labs) > 0:
            normal_times_str = " ".join([str(x) for x in sorted(normal_labs)])
            normal_res.append([pati_abrv, table_abrv_name + "." + lab_id + "+n", str(count_normal),
                               normal_times_str, "pati." + table_abrv_name])
        if len(abnormal_labs) > 0:
            abnormal_times_str = " ".join([str(x) for x in sorted(abnormal_labs)])
            abnormal_res.append([pati_abrv, table_abrv_name + "." + lab_id + "+ab", str(count_abnormal),
                                 abnormal_times_str, "pati." + table_abrv_name])


    normal_df = pd.DataFrame(normal_res, columns= file_header)
    abnormal_df = pd.DataFrame(abnormal_res, columns= file_header)
    normal_df.to_csv(os.path.join(output_dir, "Edge.csv"), index_label='NODE1', mode="a", header=False, index=False)
    abnormal_df.to_csv(os.path.join(output_dir, "Edge.csv"), index_label='NODE1', mode="a", header=False, index=False)


def write_static_event(table_events, output_dir, table_abrv_name):
    count_added = table_events.groupby(['SUBJECT_ID', 'HADM_ID', 'EVENT_ID']).size().reset_index(name="COUNT")
    count_added['NODE1'] = "pati." + count_added['SUBJECT_ID'].map(str) \
                           + "+" + count_added['HADM_ID'].map(str)
    count_added['NODE2'] = table_abrv_name + "." + count_added['EVENT_ID'].map(str)
    count_added['HOURS'] = pd.Series(["" for i in range(count_added.shape[0])],
                                     index=count_added.index)
    count_added['TYPE'] = pd.Series(["pati." + table_abrv_name for i in range(count_added.shape[0])],
                                    index=count_added.index)
    final_table = count_added[['NODE1', 'NODE2', 'COUNT', 'HOURS', 'TYPE']].drop_duplicates()
    final_table.to_csv(os.path.join(output_dir, "Edge.csv"), index_label='NODE1', mode="a", header=False, index=False)


#types = "low", "high", "phen"
def write_diagnoses(root_patient_folder, hadm_id, output_dir, type ="low", pheno_map = None):
    if type == "low":
        events = read_diagnoses(root_patient_folder)
        events = events.loc[events["HADM_ID"] == hadm_id]
        events["EVENT_ID"] = events['ICD9_CODE']
        events = events[['SUBJECT_ID', 'HADM_ID', 'EVENT_ID']]
        write_static_event(events, output_dir, "diag")


def process_episode(output_dir, root_patient_folder, episode_ind, label_type, pheno_map = None):
    try:
        episode_folder = os.path.join(root_patient_folder, "episode" + str(episode_ind + 1))
        episode = pd.read_csv(os.path.join(root_patient_folder,"episode" + str(episode_ind + 1) + ".csv"))
    except:
        print("episode csv file or folder missing for subject :" + str(root_patient_folder))
        return
    stays = read_stays(root_patient_folder)
    try:
        subject_id = stays["SUBJECT_ID"].loc[stays["ICUSTAY_ID"] == episode.Icustay.iloc[0]].iloc[0]
        hadm_id = stays["HADM_ID"].loc[stays["ICUSTAY_ID"] == episode.Icustay.iloc[0]].iloc[0]
    except:
        print("cannot match subject or admission for subject :" + str(root_patient_folder))
        return

    pati_abrv = "pati." + str(subject_id) + "+" + str(hadm_id)
    write_bio_edges(episode.Age[0], episode.Gender[0], episode.Ethnicity[0], pati_abrv, output_dir )
    write_diagnoses(root_patient_folder, hadm_id, output_dir, type = label_type, pheno_map=pheno_map)

    for table in filter(lambda x : ".csv " not in x , os.listdir(episode_folder)):
        table_abrv_name = table_name_map[(table.split(".")[0].split("_")[0])]

        if "timeseries" in table:
            with open(os.path.join(episode_folder, table)) as tsfile:
                ts_lines = tsfile.readlines()
            write_seq_event(ts_lines, output_dir, pati_abrv, table_abrv_name)
        else:
            table_events = dataframe_from_csv(os.path.join(episode_folder, table), index_col=None)
            write_static_event(table_events, output_dir, table_abrv_name)


def add_symptoms(edge_path, node_path, symptom_path, pati_map):
    symptoms = pd.read_csv(symptom_path, header=None, names=["NODE1", "NODE2", "COUNT", "TYPE"])
    symptoms = symptoms[symptoms["NODE1"].isin(pati_map)]
    symptoms['HOURS'] = pd.Series(["" for i in range(symptoms.shape[0])],
                                     index=symptoms.index)
    symptoms = symptoms[['NODE1', 'NODE2', 'COUNT', 'HOURS', 'TYPE']].drop_duplicates()
    symptoms.to_csv(edge_path, index_label='NODE1', mode="a", header=False, index=False)

    symptom_nodes = symptoms[['NODE2']]
    symptom_nodes['TYPE'] = pd.Series(["symp" for i in range(symptom_nodes.shape[0])],
                                     index=symptom_nodes.index)
    symptom_nodes.to_csv(node_path, index_label='NODE', mode="a", header=False, index=False)


def create_node_file(edge_path, node_path):
    node_map = {}
    patient_map = {}
    edges_df = pd.read_csv(edge_path, header=None, names = file_header)
    for index, row in edges_df.iterrows():
        node_map[row["NODE1"]] = row["NODE1"].split(".")[0]
        node_map[row["NODE2"]] =  row["NODE2"].split(".")[0]
        patient_map[row["NODE1"]] = True

    node_info = []
    for node in node_map:
        node_info.append([node, node_map[node]])

    node_df = pd.DataFrame(node_info, columns= ["NODE", "TYPE"])
    node_df.to_csv(node_path, index_label='NODE1', mode="a", header=False, index=False)

    return node_map, patient_map



def create_supervised_files(output_dir, edge_file_path, feature_types_keep, label_type_keep, node_map = None):
    feature_path = os.path.join(output_dir, "Features.csv")
    label_path = os.path.join(output_dir, "Labels.csv")

    open(feature_path, "w").close()
    open(label_path, "w").close()

    edges_df = pd.read_csv(edge_file_path, header=None, names=file_header)
    feature_edges_df = edges_df[edges_df["TYPE"].isin(feature_types_keep)]
    label_edges_df = edges_df[edges_df["TYPE"].isin(label_type_keep)]

    if node_map != None:
        feature_edges_df = feature_edges_df[feature_edges_df["NODE2"].isin (node_map)]
        label_edges_df = label_edges_df[label_edges_df["NODE2"].isin(node_map)]


    label_edges_df = label_edges_df[["NODE1", "NODE2"]]
    label_edges_df["LABEL"] = 1
    feature_edges_df.to_csv(feature_path, index_label='NODE1', mode="a", header=False, index=False)
    label_edges_df.to_csv(label_path, index_label='NODE1', mode="a", header=False, index=False)


def create_negative_samples(label_path, node_path, new_label_path, neg_num=100):
    nodes_df = pd.read_csv(node_path, header=None, names=["NODE", "TYPE"])
    org_label_df = pd.read_csv(label_path, header=None, names=["PATI", "DIAG", "LABEL"])
    all_diag = nodes_df["NODE"][nodes_df["TYPE"] == "diag"]

    final_res = []

    pati_diag_dict = collections.defaultdict(list)
    for index, row in org_label_df.iterrows():
        pati_diag_dict[row["PATI"]].append(row["DIAG"])

    for pati in pati_diag_dict:
        diag_list = pati_diag_dict[pati]
        white_listed_diag = list(filter(lambda a: a not in diag_list, all_diag))
        for diag in diag_list:
            final_res.append([pati, diag, 1])
        neg_samples = random.sample(white_listed_diag, neg_num - len(diag_list))
        for n in neg_samples:
            final_res.append([pati, n, 0])

    label_df = pd.DataFrame(final_res, columns=["PATI", "DIAG", "LABEL"])
    label_df.to_csv(new_label_path, index_label='PATI', mode="a", header=False, index=False)


def process_partition(args, partition):
    # preparing the edge and node pathes
    output_dir = os.path.join(args.output_path, partition)
    edge_file =  os.path.join(output_dir, "Edge.csv")
    node_file = os.path.join(output_dir, "Node.csv")
    working_dir = os.path.join(args.root_path, partition)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    open(edge_file, "w").close()
    open(node_file, "w").close()

    # iterating over patients
    patients = list(filter(str.isdigit, os.listdir(working_dir)))
    for (patient_index, patient) in enumerate(patients):
        if patient_index % 1000 == 0:
            print ("processing patient: " + str(patient_index))
        # if patient_index == 300:
        #     break
        patient_folder = os.path.join(working_dir, patient)
        for episode_ind, episode in enumerate(filter(lambda x: "episode" in x and ".csv" not in x, os.listdir(patient_folder))):
            process_episode(output_dir, patient_folder , episode_ind, args.label_type, None)

    node_map, pati_map = create_node_file(edge_file, node_file)
    add_symptoms(edge_file ,node_file, args.symptoms_file , pati_map)
    return  node_map, edge_file



def main():

    parser = argparse.ArgumentParser(description="Create data for network embedding task.")
    parser.add_argument('root_path', type=str, help="Path to root folder containing train and test sets.")
    parser.add_argument('output_path', type=str, help="Directory where the created data should be stored.")
    parser.add_argument('--label_type', '-l', type=str, default='low', help='types:low,high,pheno')
    parser.add_argument('--symptoms_file', type=str,
                        default=os.path.join(os.path.dirname(__file__), '../resources/symptoms.csv'),
                        help='CSV containing symptoms for patients.')

    parser.add_argument('--valid_suprv_types', '-t', type=str, nargs='+', help='Valid path types for supervsied model.',
                        default=['pati.bio', 'pati.labt', 'pati.proc', 'pati.symp', 'pati.micro'])
    parser.add_argument('--phenotype_definitions', '-p', type=str,
                        default=os.path.join(os.path.dirname(__file__), '../resources/hcup_ccs_2015_definitions.yaml'),
                        help='YAML file with phenotype definitions.')
    args, _ = parser.parse_known_args()


    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    print("start processing train partition")
    node_map, train_edge_file = process_partition(args, "train")
    print("making supervised files")
    create_supervised_files(os.path.join(args.output_path, "train"), train_edge_file, feature_types_keep= args.valid_suprv_types,
                            label_type_keep=["pati.diag"])
    ##############add metapath to edge file ###################
    print("start processing test partition")
    _, test_edge_file = process_partition(args, "test")
    print("making supervised files")
    create_supervised_files(os.path.join(args.output_path, "test"), test_edge_file, feature_types_keep= args.valid_suprv_types,
                            label_type_keep=["pati.diag"], node_map= node_map)
    #add negative samples
    create_negative_samples(os.path.join(args.output_path, "test","Labels.csv"),
                            os.path.join(args.output_path, "train", "Node.csv"),
                            os.path.join(args.output_path, "test", "TestLabels.csv"),
                            neg_num=100)






if __name__ == '__main__':
    main()


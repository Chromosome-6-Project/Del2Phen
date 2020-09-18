
class Patient:
    def __init__(self, ID, questionnaire):
        self.id = ID
        self.questionnaire = questionnaire

def read_genotype_data(data_file):
    with open(data_file) as infile:
        header = infile.readline()
        data = infile.readlines()
    header = header.strip("\n").split("\t")
    data = [row.strip().split("\t") for row in data]
    return header, data

def read_phenotype_data(data_file):
    with open(data_file) as infile:
        data = infile.readlines()
    data = [row.strip("\n").split("\t") for row in data]
    data = [list(x) for x in zip(*data)]
    data = [row for row in data if "".join(row)]
    data = data[1:]
    for i, row in enumerate(data):
        for j, entry in enumerate(row):
            if not entry:
                data[i][j] = None
    header = data[0]
    data = [dict(zip(header, row)) for row in data[1:]]
    return data

# def build_patients(genos, phenos):
#     id_list = genos

        
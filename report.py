class Report:

    def __init__(self):
        self.itemset = None
        self.ig = None
        self.genomes_list = []
        self.num_of_label1 = 0
        self.num_of_label2 = 0

    @staticmethod
    def print_headline(report_file_path):
        with open(report_file_path, 'w') as file:
            file.write("maxIG itemset REPORT: \n\n")

    def print_report(self, report_file_path):
        with open(report_file_path, 'a') as file:
            file.write("itemset:  %s\n" % self.itemset)
            file.write("ig:  %s\n" % self.ig)
            file.write("genomes list:  %s\n" % self.genomes_list)
            file.write("num_of_label1:  %d\n" % self.num_of_label1)
            file.write("num_of_label2:  %d\n" % self.num_of_label2)
            file.write("\n")

    def set_itemset(self, itemset):
        self.itemset = itemset

    def set_ig(self, ig):
        self.ig = ig

    def add_genome(self, genome):
        self.genomes_list.append(genome)

    def increase_num_of_label1(self):
        self.num_of_label1 += 1

    def increase_num_of_label2(self):
        self.num_of_label2 += 1
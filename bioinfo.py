import FPgrowth
import os
from math import log
import report

HABITAT_FILE = os.getcwd() + "\\bactTaxa_Habitat.txt"
TAXA_FILE = os.getcwd() + "\\taxa.txt"
COG_INFO_FILE = os.getcwd() + "\\COG_INFO_TABLE.txt"
COG_WORDS_BAC_FILE = os.getcwd() + "\\cog_words_bac.txt"
TAXA_FILE_OF_CLASSES = os.getcwd() + "\\taxa_of_classes.txt"
REPORT_FILE_PATH = os.getcwd() + "\\positive_strand_report.txt"
REPORT2_FILE_PATH = os.getcwd() + "\\negative_strand_report.txt"
LABEL1 = "Animal"
LABEL2 = "Plant"
MIN_SUP = 70 #the minimum MIN_SUP that ran!

class UtilityMethods:
    @staticmethod
    def InformationGain(patternSupport,labelSupport,patternLabelUnionSupport):
        if patternSupport == 0 or labelSupport == 0 or patternLabelUnionSupport == 0:
            return 0
        else :
            return UtilityMethods.__InformationGainFormula(patternSupport,labelSupport,patternLabelUnionSupport)

    @staticmethod
    def __InformationGainFormula(o,p,q):
        if o==0 :
            return 0

        conditionalProb_term1 = -o*q*(UtilityMethods.__Log2(q))-o*(1-q)*(UtilityMethods.__Log2(1-q))
        conditionalProb_term2 = (o*q-p)*(UtilityMethods.__Log2_wDivision((p-o*q),(1-o)))
        conditionalProb_term3 = (o*(1-q)-(1-p))*(UtilityMethods.__Log2_wDivision(((1-p)-(o*(1-q))),(1-o)))

        conditionalProb = conditionalProb_term1 + conditionalProb_term2 + conditionalProb_term3

        #use binary entropy formula to calculate entropy of the pattern
        #see http://en.wikipedia.org/wiki/Binary_entropy_function
        nonCondProb = -p*(UtilityMethods.__Log2(p))-(1-p)*(UtilityMethods.__Log2(1-p))

        # print("o: {}, p: {}, q: {}".format(o, p, q))
        #calculate and return information gain
        ig = nonCondProb - conditionalProb
        if ig < 0:
            ig = 0
        return ig

    @staticmethod
    def __Log2(x):
        ans = 0
        try:
            if x != 0:
                ans = log(x,2)
        except ValueError:
            pass
        return ans

    @staticmethod
    def __Log2_wDivision(x,y):
        ans = 0
        try:
            if x != 0 and y != 0 and (x/y) != 0:
                ans = log(x/y,2)
        except ValueError:
            pass
        return ans

def create_taxa_of_classes():
    """
    this function creates bactTaxa_Habitat file that include only the classes we chose (Animal and Plant)
    :return: 
    """
    newfile = open(TAXA_FILE_OF_CLASSES, 'w')
    with open(HABITAT_FILE) as openfileobject:
        for line in openfileobject:
            if line[line.rfind(";")+1:-1] == LABEL1 or line[line.rfind(";")+1:-1] == LABEL2:
                newfile.write(line)
    newfile.close()

def parseFile(file_path, delimiter):
    """
    this function takes a text file in which the rows are separated by newline (\n)
     and the cols are separated by delimiter given and creates a table of the original elements
    :param file_path: the path of the file we want to parse
    :param delimiter: 
    :return: table of the elements
    """
    f = open(file_path, "r")
    table = []
    text = f.read()
    count_new_lines = text.count("\n")
    count_cols = text[:text.find("\n")].count(delimiter)
    f.seek(0) # go back to the beginning of the file
    for num_line in range(count_new_lines):
        cols = []
        line = f.readline()
        for num_col in range(count_cols):
            semi_index = line.find(delimiter)
            cols.append(line[:semi_index])
            line = line[semi_index+1:]
        line = line[:line.find("\n")] # removes the \n of new line
        cols.append(line) # add the last column
        table.append(cols)
    return table

def create_table_of_chosen_habitats():
    """
    this function parses the the TAXA_FILE_OF_CLASSES to a table
    :return: chosen_habitats_table
    """
    habitats_table = parseFile(TAXA_FILE_OF_CLASSES, ";")
    chosen_habitats_table = []
    for line in habitats_table:
        chosen_habitats_table.append(line)

    return chosen_habitats_table

def create_words_bac_dict(cog_words_bac_table):
    """
    this function converts the cog_words_bac_table to cog_words_bac_dict
    the dict format is {uid: {cog_num: [words]}}
    :param cog_words_bac_table: 
    :return: cog_words_bac_dict
    """
    cog_words_bac_dict = {}
    for row_cog in cog_words_bac_table:
        if (cog_words_bac_dict.get(row_cog[3])) == None : # no such key
            cog_words_bac_dict[row_cog[3]] = {}
        cog_words_bac_dict[row_cog[3]][row_cog[0]] = row_cog[-1]
    return  cog_words_bac_dict

def create_words_bac_dict_by_strand(cog_words_bac_table, strand):
    """
    this function creates dict of words from a given strand
    :param cog_words_bac_table: 
    :param strand: 
    :return: cog_words_bac_dict_by_strand
    """
    counter = 0
    cog_words_bac_dict_by_strand = {}
    for row_cog in cog_words_bac_table:
        if row_cog[2] == strand:
            counter += 1
            if row_cog[3] not in cog_words_bac_dict_by_strand.keys() : # no such key
                cog_words_bac_dict_by_strand[row_cog[3]] = []
            cog_words_bac_dict_by_strand[row_cog[3]].append(row_cog[-1])
    return cog_words_bac_dict_by_strand

def create_dataSet(chosen_habitats_table, cog_words_bac_dict):
    """
    this function takes transactions (of cog words) from cog_words_bac_dict 
    where UID is equal to UID in chosen_habitats_table and creates a dataset of transactions
    it also creates a corresponding list of pairs (uid, label) for each transaction
    :param chosen_habitats_table: 
    :param cog_words_bac_dict: 
    :return: dataSet, labels_and_uids_list
    """
    f = open(os.getcwd() + '\\classifier.txt', 'w')
    dataSet = []
    labels_and_uids_list = [] # list of pairs (label, uid)
    for row in chosen_habitats_table:
        transaction = []
        uid = row[1]  # get the uid field
        if cog_words_bac_dict.get(uid) is not None:
            labels_and_uids_list.append((row[-1], uid)) # row[-1] informs habitat information (animal\plant)
            words_list = cog_words_bac_dict.get(uid) # the values contain the words in specific row in cog_words
            for word in words_list:
                cogs = word.split("\t")
                for cog in cogs[1:]: # the first cog is the id of the strain and the last cog will be just "\t"
                    # f.write("%s: %s" % (row[-1], cog_info_dict.get(cog)))
                    if cog not in  ('', 'X'):
                        transaction.append(cog)
            dataSet.append(transaction)
    f.close()
    print("number of transactions: ", len(dataSet))
    return dataSet, labels_and_uids_list

def count_habitats_in_dataset(labels_and_uids_list, class_label):
    """
    this function counts occurrences of 'class_label' in our dataset
    :param labels_and_uids_list: 
    :param class_label: 
    :return: count_class_label
    """
    count_class_label = 0
    for pair in labels_and_uids_list:
        if pair[0] == class_label:
            count_class_label+=1
    return count_class_label

def pre_processing_IG(dataset, labels_and_uids_list, itemset, label):
    """
    this function creates the variables needed for IG computation
    o = the frequency of the 'itemset' in dataset
    p = the frequency of 'label' in dataset
    q = the frequency of 'label' in the transactions containing 'itemset'
    :param dataset: 
    :param labels_and_uids_list: 
    :param itemset: 
    :param label: 
    :return: 
    """
    o = 0
    p = 0
    q = 0
    counter_subset_occurrences = 0
    itemset_label_union_occurrences = 0
    for i, trans in enumerate(dataset):
        set_trans = set(trans)
        if itemset.issubset(set_trans):
            counter_subset_occurrences += 1
            if labels_and_uids_list[i][0] == label:
                itemset_label_union_occurrences += 1
    num_of_transactions = len(dataset)

    if num_of_transactions != 0:
        number_of_label = count_habitats_in_dataset(labels_and_uids_list, label)
        p = number_of_label / num_of_transactions
        o = counter_subset_occurrences / num_of_transactions
    if counter_subset_occurrences != 0:
        q = itemset_label_union_occurrences / counter_subset_occurrences
    return o, p, q

def get_max_IG(dataset, labels_and_uids_list, freqItemList):
    """
    this function finds the itemset from 'freqItemList' with the maximum IG, and returns both the maxIG and the itemst
    with this ig.
    :param dataset: 
    :param labels_and_uids_list: 
    :param freqItemList: 
    :return: max_IG, itemset_with_max_IG
    """
    max_IG = 0
    itemset_with_max_IG = None
    for itemset in freqItemList:
        print("itemset: ", itemset)
        itemset_support, label_support, itemset_label_union_support = pre_processing_IG(dataset, labels_and_uids_list, itemset, LABEL1)
        ig = UtilityMethods.InformationGain(itemset_support, label_support, itemset_label_union_support)
        print("ig: ", ig)
        if ig > max_IG:
            max_IG = ig
            itemset_with_max_IG = itemset
    return max_IG, itemset_with_max_IG

def remove_trans_with_itemset(itemset, dataset, labels_and_uids_list):
    """
    this function removes transactions that contain 'itemset' and their label and uid from the
    corresponding 'labels_and_uids_list'
    :param itemset: 
    :param dataset: 
    :param labels_and_uids_list: 
    :return: dataset, labels_and_uids_list
    """
    i = len(dataset)-1
    for trans in reversed(dataset):
        set_trans = set(trans)
        if itemset.issubset(set_trans):
            dataset.remove(trans)
            labels_and_uids_list.remove(labels_and_uids_list[i])
        i -= 1
    return dataset, labels_and_uids_list

def print_report(itemset, ig, dataset, labels_and_uids_list, label):
    """
    this function write details about the given 'itemset' to the report file
    :param itemset: 
    :param ig: 
    :param dataset: 
    :param labels_and_uids_list: 
    :param label: 
    :return: 
    """
    rep = report.Report()
    rep.set_itemset(itemset)
    rep.set_ig(ig)
    for i, trans in enumerate(dataset):
        set_trans = set(trans)
        if itemset.issubset(set_trans):
            rep.add_genome(labels_and_uids_list[i][1])
            if labels_and_uids_list[i][0] == label:
                rep.increase_num_of_label1()
            else:
                rep.increase_num_of_label2()

    rep.print_report(REPORT_FILE_PATH)

if __name__ == '__main__':
    freqItemList = []
    # parsing HABITAT_FILE to 'animals_and_plants_table'
    animals_and_plants_table = create_table_of_chosen_habitats()
    # Parsing COG_WORDS_BAC_FILE to 'cog_words_bac_dict'
    cog_words_bac_table = parseFile(COG_WORDS_BAC_FILE, "#")

    positive_strand_dict = create_words_bac_dict_by_strand(cog_words_bac_table, '1')
    negative_strand_dict = create_words_bac_dict_by_strand(cog_words_bac_table, '-1')
    cog_words_bac_dict = create_words_bac_dict(cog_words_bac_table)
    strands_dicts_list = [positive_strand_dict, negative_strand_dict]

    # creating the dataset of each strand
    for strands_dict in strands_dicts_list:
        dataset, labels_and_uids_list = create_dataSet(animals_and_plants_table, strands_dict)
        report.Report.print_headline(REPORT_FILE_PATH)

        # main loop, will run on the dataset and the updated dataset until no transaction is left OR the FPtree is empty
        while(len(dataset) > 0):
            # Creating the FPtree of dataset
            initSet = FPgrowth.createInitSet(dataset)
            myFPtree, header_table = FPgrowth.createTree(initSet, MIN_SUP)
            print("the tree was created!")

            if myFPtree is not None:
                # finding freqItemList
                FPgrowth.mineTree(myFPtree, header_table, MIN_SUP, set([]), freqItemList)
                max_IG, itemset_with_max_IG = get_max_IG(dataset, labels_and_uids_list, freqItemList)
                print("itemset: {}, max IG: {}".format(itemset_with_max_IG, max_IG))
                if itemset_with_max_IG is None:
                    break
                print_report(itemset_with_max_IG, max_IG, dataset, labels_and_uids_list, LABEL1)
                # updating the dataset
                dataset, labels_and_uids_list = remove_trans_with_itemset(itemset_with_max_IG, dataset, labels_and_uids_list)
            else: # no transaction met the MIN_SUP
                break
            freqItemList = []

        REPORT_FILE_PATH = REPORT2_FILE_PATH
        print("End of strand")
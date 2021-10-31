# bacterias-pattern-distinguish
This project was done within the framework of the course "Topics in Bioinformatic" by Prof. Ziv-ukelson Michal, graded 100


**Summary of the steps in the project:**

We were given a dataset of bacteria that were annotated with their environmental habitat. Then we were asked to select two of the habitats for our project, we chose the bacteria that their habitats are "Animals" or "Plants".

For each bacterium in our data we were given a set of COGS (Clusters of Orthologous Groups) words, each set was called a transaction.

**The main purpose was to learn from the data and find subsets (itemsets) from the transaction that can, with high probably, distinguish between animal bacteria and plant bacteria.**

We started by filtering and parsing the given data to be able to create an FP tree structure for all the transactions of the bacteria in our filtered data. 
Then we used the code for the FP-Growth algorithm to generate a program to solve the following problems:

Problem 1: 
Extract all frequent itemsets from the transactions (using the FP-Growth algorithm) with support >= min_sup and compute IG (Information Gain) for each frequent itemset. 
Sort the frequent itemsets by decreasing IG score value and identify the itemset that has the highest IG score. 
If the search is too slow to complete in reasonable time, use a very large min_sup value to constrain it.

Problem 2: 
Run the program you generated in problem 1 and find highest scoring Discriminative Itemset 𝑫𝜶 in D. 
Update the transactional data set D by removing any transaction that contains 𝑫𝜶 and update the FP tree accordingly. 
Repeat until no transaction is left or no itemset passes min_sup. 
The output of this step consists of all itemsets that were identified in one of the iterations as MAX_IG. 

Because the data we were given was very large we did not receive a good and informative result, therefore, on the second part of the project we were asked to repeat the experiment, but this time “invent” our own additional constraints (either biological or computational) that will either improve the biological reasoning or speed up the search while maintaining its accuracy.

**To learn about the additional constraints that we suggested and our biological conclusions from the algorithm results see "Report Part C".**

Note: To run the program, there is a need to unzip the cog_words_bac file, by the command 'gunzip data/cog_words_bac.txt.gz'

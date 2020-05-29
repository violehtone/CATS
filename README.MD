# CATS (Classification Assessment of Tumor Subtypes) assignment
Assignment for the course Bioinformatics for translational medicine
M.Sc. in Bioinformatics & Systems biology @ VU / UvA


### Introduction
In this group assignment, you need to train a classifier that is able to predict subtypes of breast cancer
tumors based on array CGH (aCGH) data.

### Background
Breast cancer is a heterogeneous disease and classification of breast cancer tumors in their molecular
subtypes has important implications for the treatment and prognosis. Three receptors play a pivotal role in
these subtypes: the Estrogen Receptor (ER), Progesterone Receptor (PR) and Human Epidermal growth
factor Receptor 2 (HER2). After removal of the tumor, the pathology department of the hospital tests these
samples for presence of ER, PR and HER2. The three main subtypes in breast cancer on which the
treatment decision will be based are:

- HER2 positive: HER2+
- Hormone receptor positive (HR+): ER+ and/or PR+, and HER2-
- Triple negative (TN): ER-, PR- and HER2-

Each of the three subtypes reacts differently to different types of treatment.

### Data
We provide you with a dataset with 100 breast cancer samples from the three subtypes. These samples
are analyzed on a high-resolution array CGH platform with 244,000 probes per array that measures the
quantity of chromosomal DNA. The pre-processing of the data has been done for you.
For each of these regions and each sample we give the call whether that region is a gain, an amplification,
a loss or normal (-1 for loss, 0 for normal, 1 for gain, and 2 for amplification of the DNA). Two files are
provided, one containing the preprocessed aCGH data of the cancer samples (Train_call.txt) and the other
containing the associated clinical outcome of these samples — the subtypes (Train_clinical.txt).

### Building your classifier
Your task is to get a well-trained classifier for predicting the three breast cancer subtypes.
To achieve this, you might need some of the following steps.
1. Data purification, transformation, if necessary.
2. Feature selection, if necessary.
3. Choose machine learning methods (classifiers).
4. Train and validate the classifiers.
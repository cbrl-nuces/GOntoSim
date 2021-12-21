# GOntoSim: A Semantic Similarity Measure based on LCA and Common Descendants

Web server is available at http://www.cbrlab.org/GOntoSim.html.

## Dependencies
- Goatools:
Make sure your Python version >= 2.7, install the latest stable version via PyPI:

```
pip install goatools
```



## Data
Annotations for these annotations have already been downloaded and filtered. These are saved as named tuples required in the code. 

- 'ALL_EC_GOTERMS_IEA_BP.py'
- 'ALL_EC_GOTERMS_IEA_MF.py'
- 'ALL_EC_GOTERMS_NonIEA_BP.py'
- 'ALL_EC_GOTERMS_NonIEA_MF.py'

The Gene Ontology File (go-basic.obo) used in the experiments is also provided in the same folder. 

The association files are for the implementation of Resnik's and Lin's measures.
- 'association_IEA_MF.txt'
- 'association_nonIEA_MF.txt'

## iPython Notebooks

- 'Semantic Similarity-code.ipynb' has the complete implementation. 

- 'GOntoSim-Example Usage.ipynb' has a few examples for calculating similarities between 2 given GO-Terms.

- Download GO annotations And Convert to Named Tuples.ipynb contains the code to download GO Annotations using the QuickGO API, filter the annotations as required, 
and contains the code to convert the annotations to the named tuples required in the code.

## To replicate the complete set of Experiments

Any one of the following measures can be used to calculate the similarity:
'gontosim',	'baseline', 'lca', 'baselineDesc', 'wang', 'gogo' 
('resnik' and 'lin' can be calculated for MF only)

The arguments required are the similarity measure, the GO Aspect (MF or BP), Evidence Code (IEA or NONIEA), number of samples 
```
python GOntoSim.py measure GO_Aspect Evidence_Code Number_of_Samples

```

Run GOntoSim.py with the following Commands:

### Experiment 1: 
	
This experiment uses the Molecular Function GO term annotations (NONIEA) for the Enzymes.  
	
```
python GOntoSim.py gontosim MF NONIEA 150

```
### Experiment 2: 
	
This experiment uses the Biological Process GO term annotations (NONIEA) for the Enzymes.  
	
```
python GOntoSim.py gontosim BP NONIEA 150

```
### Experiment 3: 
	
This experiment uses the Molecular Function GO term annotations (IEA) for the Enzymes.  
	
```
python GOntoSim.py gontosim MF IEA 500

```






Contact hammad.naveed@nu.edu.pk or amna.kamran@nu.edu.pk

Feel free to contact us in case of any confusions.

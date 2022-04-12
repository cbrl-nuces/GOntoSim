import GOntoSim
from goatools.obo_parser import GODag
from goatools.base import get_godag
go = get_godag("go-basic.obo", optional_attrs={'relationship'})


def main():
	gene1 = ['GO:0004022','GO:0004024', 'GO:0004174', 'GO:0046872','GO:0008270','GO:0004023', 'GO:0016491']
	gene2 = ['GO:0009055','GO:0005515','GO:0046872','GO:0008270','GO:0020037']
	
	# getting a list of unique GO Terms
	unique_goterms = list(set().union(gene1, gene2))
	
	method = 'GOntoSim'	
	# Setting the variable S_values
	if method == 'GOntoSim':
		S_values = [(x,GOntoSim.Semantic_Value(x, go, 'Baseline_LCA_avg')) for x in unique_goterms ]
	else:
		print(method)
		S_values = [(x,GOntoSim.Semantic_Value(x, go, method)) for x in unique_goterms ]
	
	S_values = dict(S_values)

	print(method)
	print("Similarity: ", GOntoSim.Similarity_of_Set_of_GOTerms(gene1, gene2, method, S_values))

if __name__ == '__main__':
	main()

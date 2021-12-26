## import Libraries

from goatools.semantic import common_parent_go_ids, min_branch_length, semantic_distance, semantic_similarity
from goatools.obo_parser import GODag
from goatools.base import get_godag
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim
from goatools.semantic import TermCounts, get_info_content
from goatools.associations import read_associations
from goatools.semantic import deepest_common_ancestor

import numpy as np
import pandas as pd
from collections import namedtuple
from collections import defaultdict

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn import metrics
from scipy.cluster.hierarchy import dendrogram 
import scipy.cluster.hierarchy as shc
from scipy.spatial.distance import squareform

import sys
import random

# import annotations
from ALL_EC_GOTERMS_IEA_MF import *
from ALL_EC_GOTERMS_NonIEA_MF import *
from ALL_EC_GOTERMS_IEA_BP import *
from ALL_EC_GOTERMS_NonIEA_BP import *

# define named tuple
Gene = namedtuple('Gene', 'GeneName GOAnnotations')
# load Gene Ontology
go = get_godag("go-basic.obo", optional_attrs={'relationship'})
# Define Relationship Semantic Contribution 
is_a = 0.8
part_of = 0.6
regulates = 0.7
negatively_regulates = 0.7
positively_regulates = 0.7
reg = 'reg0.7'
c = 0.67    #for GOGO
## helper functions
def all_common_parent_go_ids(goids, godag):
	'''
		This function finds the common ancestors in the GO
		tree of the list of goids in the input.
	'''
	# Find candidates from first
	rec = godag[goids[0]]
	candidates = rec.get_all_upper()
	candidates.update({goids[0]})

	# Find intersection with second to nth goid
	for goid in goids[1:]:
		rec = godag[goid]
		parents = rec.get_all_upper()
		parents.update({goid})

		# Find the intersection with the candidates, and update.
		candidates.intersection_update(parents)
	return candidates
def lowest_common_ancestor(goterms, godag):
	'''
		This function gets the nearest common ancestor
		using the above function.
		Only returns single most specific - assumes unique exists.
	'''
	# Take the element at maximum depth.
	return max(all_common_parent_go_ids(goterms, godag), key=lambda t: godag[t].depth)
def all_paths_to_top(term, godag):
	# inputs: term_id and Go dag with 'relationship' as optional attributes
		""" Returns all possible paths to the root node"""
		if term not in godag:
			sys.stderr.write("Term %s not found!\n" % term)
			return
		def _all_paths_to_top_recursive(rec):
			if rec.level == 0:
				return [[rec]]
			paths = []
			parents = rec.get_goterms_upper()
			for parent in parents:
				top_paths = _all_paths_to_top_recursive(parent)
				for top_path in top_paths:
					top_path.append(rec)
					paths.append(top_path)
			return paths

		go_term = godag[term]
		return _all_paths_to_top_recursive(go_term)
    
    
def all_paths_to_top_wang(term, godag, optional_relationships):
	# inputs: term_id and Go dag with 'relationship' as optional attributes
		""" Returns all possible paths to the root node"""
		if term not in godag:
			sys.stderr.write("Term %s not found!\n" % term)
			return
		def _all_paths_to_top_recursive(rec):
			if rec.level == 0:
				return [[rec]]
			paths = []
			parents = rec.get_goterms_upper_rels(optional_relationships)
			#parents = rec.get_goterms_upper()
			for parent in parents:
				#parent = go[parent1]
				top_paths = _all_paths_to_top_recursive(parent)
				for top_path in top_paths:
					top_path.append(rec)
					paths.append(top_path)
			return paths
		go_term = godag[term]
		return _all_paths_to_top_recursive(go_term)
def Semantic_Value(go_id, go, method):
	''' input: goterm_id 
		returns all of the weighted (all relatinships in go-basic) paths to root 
		#relationship types are global variables with appropriate weights
	'''
	if method == 'wang':
		# calculates all paths to top (with all relationships)
		optional_relationships = {'part_of'}
		all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
		S_values = list()
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			#print("index = " + str(index))
			path.reverse()
			for idx, term in enumerate(path):
				# Semantic Value of the term itself is always 1 
				if idx == 0:	#which is on idx = 0 
					S_values[index].append((go_id, 1))
				if idx < len(path)-1:
					if term.relationship != {}:
						#print(term.relationship.keys())
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
							else: 
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						else: 
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
					else: 
						S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
		return final_values(S_values,'max')
	elif method == 'GOGO':
		optional_relationships = {'part_of'}
		all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
		S_values = list()
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			#print("index = " + str(index))
			path.reverse()
			for idx, term in enumerate(path):
				# Semantic Value of the term itself is always 1 
				if idx == 0:	#which is on idx = 0 
					S_values[index].append((go_id, 1))
				if idx < len(path)-1:
					if term.relationship != {}:
						#print(term.relationship.keys())
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
								weight = (1 / (c + len(go[path[idx+1].item_id].children))) + part_of
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight))
							else: 
								weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
						else: 
							weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
					else: 
						weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
						S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
		return final_values(S_values,'max')
	#Baseline Measure - Almost Same as Semantic Value; only Difference is the weight of root as ancestor = 0 and different realtionships
	else:
		all_all_paths = all_paths_to_top(go_id, go)
		S_values = list()
		#print(go_id)
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			#print("index = " + str(index))
			path.reverse()
			for idx, term in enumerate(path):
				# Semantic Value of the term itself is always 1 
				if idx == 0:	#which is on idx = 0 
					S_values[index].append((go_id, 1))
				if idx < len(path)-1 and path[idx+1].item_id != 'GO:0003674' and path[idx+1].item_id != 'GO:0005575' and path[idx+1].item_id != 'GO:0008150':
					if term.relationship != {}: 
						#print(term.relationship.keys())
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
							else: 
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'regulates' in term.relationship:
							if path[idx+1] in term.relationship['regulates']:
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * regulates))
							else: 
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'negatively_regulates' in term.relationship:
							if path[idx+1] in term.relationship['negatively_regulates']:
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * negatively_regulates))
							else: 
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'positively_regulates' in term.relationship:
							if path[idx+1] in term.relationship['positively_regulates']:
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * positively_regulates))
							else: 
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
						else: 
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					else: 
						S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
				if (term.item_id == 'GO:0003674' or term.item_id == 'GO:0005575' or term.item_id == 'GO:0008150'):
						S_values[index].append((term.item_id, 0))
		if method == 'Baseline':
			return final_values(S_values, 'max')
		if method == 'Baseline_LCA' :
			svaluesfinal = final_values(S_values, 'max')
			# Replace the max or min s-values in all paths to assign only one value to each node, consistent in all paths
			S_values_Modified = list()
			for index, path in enumerate(S_values):
				S_values_Modified.append([])
				for idx, term1 in enumerate(path):
			#		print(list(value for (term, value) in svaluesfinal if term == term1[0]))
					ind = [value for (term, value) in svaluesfinal if term == term1[0]]
					S_values_Modified[index].append((path[idx][0],ind[0]))
			SumOfNodesOnEachPath = list()
			for index, path in enumerate(S_values_Modified):
			#	print("Modified index and path: ", index, path)
				SumOfNodesOnEachPath.append(sum(x[1] for x in path))
			maxPath = SumOfNodesOnEachPath.index(max(SumOfNodesOnEachPath))
			return S_values_Modified[maxPath], SumOfNodesOnEachPath[maxPath]

def final_values(S_values, isMax):
	''' helper function to assign the max of the weights assigned to each term'''
	#S_values = sorted(S_values, key=lambda x: x[0])
	unique_terms_s_values = []
	for path in S_values:
		for term in path:
			unique_terms_s_values.append(term)
	unique_terms_s_values = sorted(unique_terms_s_values, key=lambda x: x[0])
	_s_values = {}
	for y, x in unique_terms_s_values: 
		if y in _s_values: 
			_s_values[y].append((y,x)) 
		else: 
			_s_values[y] = [(y, x)]
	final_s_values = []
	if(isMax == 'max'):
		for node in _s_values:
			final_s_values.append(max(_s_values[node]))
	elif (isMax == 'min'):
		for node in _s_values:
			final_s_values.append(min(_s_values[node]))
	return final_s_values
def intersection(lst1, lst2): 
	'''Helper Function to find intersecting terms from the two input lists of (term, s_Value)'''
	da = {v:k for v, k in lst1}
	db = {v:k for v, k in lst2} 
	return [(da[k],db[k]) for k in da.keys() & db.keys()]
## Downward Graph 				  
def common_children_go_ids(goids, godag):
	'''
		This function finds the common children in the GO
		tree of the list of goids in the input.
	'''
	# Find candidates from first
	rec = godag[goids[0]]
	candidates = rec.get_all_lower()
	candidates.update({goids[0]})
	# Find intersection with second to nth goid
	for goid in goids[1:]:
		rec = godag[goid]
		children = rec.get_all_lower()
		children.update({goid})
		# Find the intersection with the candidates, and update.
		candidates.intersection_update(children)
	return candidates
def highest_common_descendant(goterms, godag):
	'''
		This function gets the nearest common descendant
		using the above function.
		Only returns single most specific - assumes unique exists.
	'''
	# Take the element at minimum depth.
	common_children = common_children_go_ids(goterms, godag)
	if len(common_children) != 0:
		# take depth attribute instead of depth to accomodate all relationships
		return min(common_children, key=lambda t: godag[t].depth)
	else:
		return 0
def all_paths_to_bottom(term, godag, x):
		# inputs: term_id and Go dag with 'relationship' as optional attributes
	""" Returns all possible paths to the root node"""
	if term not in godag:
		sys.stderr.write("Term %s not found!\n" % term)
		return

	def _all_paths_to_bottom_recursive(rec):

		if rec.depth == godag[term].depth+x:
			return [[rec]]
		else:
			paths = []
			children = rec.get_goterms_lower()
			for child in children:
				bottom_paths = _all_paths_to_bottom_recursive(child)
				for bottom_path in bottom_paths:
					bottom_path.append(rec)
					paths.append(bottom_path)
			return paths
	go_term = godag[term]
	return _all_paths_to_bottom_recursive(go_term)
def Downward_Semantic_Value(go_id, go, x):
	''' input: goterm_id 
		returns all of the weighted nodes in path to the Goterm from the children at the xth level below
		#relationship types are global variables with appropriate weights
	'''
	all_all_paths = all_paths_to_bottom(go_id, go, x)
	#print(len(all_all_paths))
	S_values = list()
	for index, path in enumerate(all_all_paths):
		S_values.append([])
		#print("index = " + str(index))
		path.reverse()
		for idx, term in enumerate(path):
			#print (term.item_id)
			if idx == 0:	#which is on idx = 0 
				S_values[index].append((go_id, 1))
			if idx < len(path)-1:
				if term.relationship != {}: 
					#print(term.relationship.keys())
					if 'part_of' in term.relationship:
						if path[idx+1] in term.relationship['part_of']:
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
						else: 
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'regulates' in term.relationship:
						if path[idx+1] in term.relationship['regulates']:
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * regulates))
						else: 
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'negatively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['negatively_regulates']:
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * negatively_regulates))
						else: 
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'positively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['positively_regulates']:
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * positively_regulates))
						else: 
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
				else: 
					S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
	return final_values(S_values, 'max')
### Calculating Similarity 
def Similarity_of_Two_GOTerms(go_id1, go_id2, go, method):
	if method == 'Baseline_LCA':
		lca = lowest_common_ancestor((go_id1, go_id2), go)
		#print(lca)
		sim_lca = Semantic_Value(lca, go, method)
		sim1 = Semantic_Value(go_id1, go, method)
		sim2 = Semantic_Value(go_id2, go, method)
		# sim1[1] and sim2[1] are the sums of all the nodes on the path with the max s-values. 
		sum_sim1_sim2 = (sim1[1] + sim2[1])
		return ((sim_lca[1]*2)/sum_sim1_sim2)
	
	elif method == 'GOntoSim':
		hcd = highest_common_descendant((go_id1, go_id2), go)
		if hcd != 0:
			hcd_depth = go[hcd].depth
			go1_depth = go[go_id1].depth
			go2_depth = go[go_id2].depth
			x = hcd_depth - go1_depth
			y = hcd_depth - go2_depth
			sv_a = Downward_Semantic_Value(go_id1, go, x)
			sv_b = Downward_Semantic_Value(go_id2, go, y)
			intersecting_terms = intersection(sv_a,sv_b)
			numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
			sim_down = (numerator/denominator)
		else:
			sim_down = 0
		sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline_LCA')
		sim = (sim_down*0.5) + (sim_upper*0.5)
		#sim = (sim_down*0.3) + (sim_upper*0.7)
		return sim
		
	elif method == 'wang' or method == 'Baseline' or method == 'GOGO':
		sv_a = Semantic_Value(go_id1, go, method)
		sv_b = Semantic_Value(go_id2, go,  method)
		intersecting_terms = intersection(sv_a,sv_b)
		numerator = sum([x for t in intersecting_terms for x in t])
		denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
		Similarity = (numerator/denominator)
		return Similarity
	elif method == 'Baseline_Desc':
		hcd = highest_common_descendant((go_id1, go_id2), go)
		if hcd != 0:
			hcd_depth = go[hcd].depth
			go1_depth = go[go_id1].depth
			go2_depth = go[go_id2].depth
			x = hcd_depth - go1_depth
			y = hcd_depth - go2_depth
			sv_a = Downward_Semantic_Value(go_id1, go, x)
			sv_b = Downward_Semantic_Value(go_id2, go, y)
			intersecting_terms = intersection(sv_a,sv_b)
			numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
			sim_down = (numerator/denominator)
		else:
			sim_down = 0
		sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline')
		sim = (sim_down*0.5) + (sim_upper*0.5)
		#sim = (sim_down*0.3) + (sim_upper*0.7)
		return sim
	elif method == 'Baseline_Desc_only':
		hcd = highest_common_descendant((go_id1, go_id2), go)
		if hcd != 0:
			hcd_depth = go[hcd].depth
			go1_depth = go[go_id1].depth
			go2_depth = go[go_id2].depth
			x = hcd_depth - go1_depth
			y = hcd_depth - go2_depth
			sv_a = Downward_Semantic_Value(go_id1, go, x)
			sv_b = Downward_Semantic_Value(go_id2, go, y)
			intersecting_terms = intersection(sv_a,sv_b)
			numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
			sim_down = (numerator/denominator)
		else:
			sim_down = 0
		return sim_down
def Similarity_of_Set_of_GOTerms(set1, set2, method):
	Sim1 = []
	Sim2 = []
	for idx, goterm in enumerate(set1):
		# print ("=========", goterm)
		Sim1.append([])
		for goid in set2:
			#print (goterm , goid)
			Sim1[idx].append((goterm, goid,(Similarity_of_Two_GOTerms(goterm, goid, go, method))))
	for idx, goterm in enumerate(set2):
		Sim2.append([])
		for goid in set1:
			#print (goterm , goid)
			Sim2[idx].append((goterm, goid,(Similarity_of_Two_GOTerms(goterm, goid, go, method))))
	sem1 = []
	sem2 = []

	for index, goterm in enumerate(Sim1):	
		sem1.append((max(Sim1[index], key=lambda x: x[2])))
	for index, goterm in enumerate(Sim2):	
		sem2.append((max(Sim2[index], key=lambda x: x[2])))

	similarity = (sum(x[2] for x in sem1)+ sum(x[2] for x in sem2) )/(len(set1) + len(set2))
	
	return round(similarity, 3)
## Evaluation: Clustering, AMI/ARI Scores

def Similarity_Matrix(genes, method, S_values):
	sim_matrix = []
	for idx,gene in enumerate(genes):
		print(gene)
		sim_matrix.append([(lambda x: Similarity_of_Set_of_GOTerms(x[1],gene[1], method, S_values))(x) for x in genes])
	return sim_matrix

def Agglomerative_Clustering(pathway, Genes, n_clusters, method, S_values):
	# Similarity Matrix
	data = Similarity_Matrix(Genes, method, S_values)
	length = len(data)
	data1 = pd.DataFrame(data = data, index=[x[0] for x in Genes],columns=[x[0]for x in Genes])
	#print('similarity matrisx: ')
	#print(data1)
	#writeSim = pathway + "_SimilarityMatrix.csv"
	#data1.to_csv(writeSim)
	data_matrix = []
	for row in data:
		data_matrix.append([(lambda x: 1-x)(x) for x in row])
	my_data = pd.DataFrame(data = data_matrix, index=[x[0] for x in Genes],columns=[x[0]for x in Genes])
	# print('distance matrix: ')
	# print( my_data)
	# writeDist = pathway + "_DistanceMatrix.csv"
	# my_data.to_csv(writeDist)	
	return Agglomerative(my_data,Genes, pathway, n_clusters)

def Agglomerative( data, Genes, pathway, n_clusters):
	model = AgglomerativeClustering(n_clusters, affinity='precomputed', linkage='complete').fit_predict(data)
	GeneNames = [gene.GeneName for gene in Genes]
	return model.tolist()

def ARI_Score(label_pred, labels_true):
	return adjusted_rand_score(label_pred, labels_true)

def AMI_Score(label_pred, labels_true):
	return adjusted_mutual_info_score(label_pred, labels_true, average_method='arithmetic')

def cont_matrix (y_true, y_pred):
	contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
	# return purity
	return(contingency_matrix)

def purity_score(y_true, y_pred):
	# compute contingency matrix (also called confusion matrix)
	contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
	#print(contingency_matrix)
	return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix) 

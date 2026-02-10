# PART 3: Exploring the gene ontology---------------

# Install the required packages
# !pip install goatools pygraphviz matplotlib 

# Load the packages
from goatools import obo_parser
import matplotlib.pyplot as plt

# Download the basic version of the gene ontology file 
# Check https://geneontology.org/docs/download-ontology/
# for more information

# Computer readable ontology file is located 
# at https://purl.obolibrary.org/obo/go/go-basic.obo
# use wget to download to your working directory


go_dict = obo_parser.GODag("data/go-basic.obo", load_obsolete = False)
# The file is loaded as a dictionary object; {key:value} pairs

# Hint, len(dictionary)
total_terms = len(go_dict)
print(f"Total number of functional active terms: {total_terms}\n")


go_dict_obsolete = obo_parser.GODag("data/go-basic.obo", load_obsolete = True)

total_terms_obsolete = len(go_dict_obsolete)
print(f"Total number of functional active + obsolete terms: {total_terms_obsolete}\n")
print(f"Total number of functional obsolete terms: {total_terms_obsolete - total_terms}\n")


with open("Figures/tot_number_of_terms.txt", "w") as f:
	f.write(f"Total number of functional active terms: {total_terms}\n")
	f.write(f"Total number of functional active + obsolete terms: {total_terms_obsolete}\n")
	f.write(f"Total number of functional obsolete terms: {total_terms_obsolete - total_terms}\n")

term_of_interest = "GO:0042710"

# textual name of functional term
go_dict[term_of_interest].name
print(f"Name of term {term_of_interest}: {go_dict[term_of_interest].name}\n")

# Depth of term
go_dict[term_of_interest].depth
print(f"Depth of term {term_of_interest}: {go_dict[term_of_interest].depth}\n")

# subcategory
go_dict[term_of_interest].namespace
print(f"Namespace of term {term_of_interest}: {go_dict[term_of_interest].namespace}\n")

with open("Figures/Term_statistics.txt", "w") as f:
	f.write(f"Name of term {term_of_interest}: {go_dict[term_of_interest].name}\n")
	f.write(f"Depth of term {term_of_interest}: {go_dict[term_of_interest].depth}\n")
	f.write(f"Namespace of term {term_of_interest}: {go_dict[term_of_interest].namespace}\n")

# draw lineage, Hint: go_dict.draw_lineage(...)
plt.figure(figsize=(12, 10))
go_dict.draw_lineage([go_dict[term_of_interest]])
plt.show()

# count #terms containing a given text
# initialize counter variable
ligase_count = 0

for term in go_dict.values():
  	if "ligase" in term.name.lower():
	  	ligase_count += 1
		  
print(f"Number of functional terms containing 'ligase': {ligase_count}")

with open("Figures/Ligase.txt", "w") as f:
	f.write(f"Number of functional terms containing 'ligase': {ligase_count}\n")

# initialize list to store the depths of all functional terms
# in the ontology
depths = []
for term in go_dict.values():
	depths.append(term.depth) 

depth_counts = {}

# 2Creating a dictionary to store depth counts of each length
for d in depths:
    if d in depth_counts:
        depth_counts[d] += 1 
    else:
        depth_counts[d] = 1


# axes labels, title, etc
plt.figure(figsize=(10, 6))
plt.bar(depth_counts.keys(), depth_counts.values(), color='skyblue', edgecolor='black')

plt.title('Ontology Width: Number of Terms per Depth')
plt.xlabel('Depth')
plt.ylabel('Number of Terms')

plt.savefig('Figures\Ontology_width_chart.png', dpi=300, bbox_inches='tight')
plt.show()








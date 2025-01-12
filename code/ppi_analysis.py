import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
from io import StringIO

# Step 1: Load the DEGs CSV
deg_file = '/content/drive/MyDrive/Mouse_AD/deg_5xFAD.csv'  # Path to your CSV file
deg_df = pd.read_csv(deg_file)

# Step 2: Extract the gene list (ENSEMBL IDs)
gene_list = deg_df['gene'].tolist()  # Use the `gene` column (ENSEMBL IDs)

# Step 3: Function to fetch PPI data from STRING
def fetch_string_ppi(gene_list, species='10090'):  # 10090 is the taxon ID for mouse
    base_url = "https://string-db.org/api/tsv/network"
    genes = "\n".join(gene_list)  # Join gene IDs into a newline-separated string
    params = {'identifiers': genes, 'species': species, 'limit': 1}
    response = requests.post(base_url, data=params)
    if response.status_code == 200:
        return response.text  # Return raw PPI data
    else:
        print("Failed to fetch data from STRING:", response.status_code)
        return None

# Step 4: Fetch PPI network data
ppi_data = fetch_string_ppi(gene_list, species='10090')
if ppi_data is None:
    raise ValueError("No PPI data retrieved. Check your input or STRING database connection.")

# Step 5: Parse the PPI data into a DataFrame
ppi_df = pd.read_csv(StringIO(ppi_data), sep="\t")

# Step 6: Filter interactions based on score > 0.7
ppi_df_filtered = ppi_df[ppi_df['score'] > 0.7]

# Step 7: Create the PPI network
# Create the PPI network
G = nx.Graph()

# Add edges to the graph using 'preferredName_A', 'preferredName_B', and 'score'
for _, row in ppi_df_filtered.iterrows():
    G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])

# Step 8: Visualize the network
plt.figure(figsize=(12, 12))
nx.draw_networkx(G, node_size=50, with_labels=False, font_size=10, width=1, alpha=0.7)
plt.title('PPI Network for DEGs in 5xFAD')
plt.show()

# print # of nodes, edges, and network density 
print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")
print(f"Network density: {nx.density(G)}")

# Find connected components
components = list(nx.connected_components(G))
print(f"Number of clusters: {len(components)}")

# Get the largest cluster
largest_component = max(components, key=len)
subgraph = G.subgraph(largest_component)

# Visualize the largest cluster
plt.figure(figsize=(10, 10))
nx.draw_networkx(subgraph, node_size=50, with_labels=False)
plt.title("Largest Cluster in PPI Network")
plt.show()

# Find top 10 nodes by degree centrality 
degree_centrality = nx.degree_centrality(G)
top_degree_nodes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:10]
print("Top 10 nodes by degree centrality:", top_degree_nodes)

# Betweenness centrality: important nodes bridging clusters
betweenness_centrality = nx.betweenness_centrality(G)

# Top nodes by betweenness
top_betweenness_nodes = sorted(betweenness_centrality.items(), key=lambda x: x[1], reverse=True)[:10]
print("Top 10 nodes by betweenness centrality:", top_betweenness_nodes)

# Calculate the clustering coefficient for each node
clustering_coefficients = nx.clustering(G)

# Sort the nodes by clustering coefficient in descending order and get the top 10
top_10_clustering_nodes = sorted(clustering_coefficients.items(), key=lambda item: item[1], reverse=True)[:10]

# Print the top 10 nodes with their clustering coefficients
print("Top 10 nodes with highest clustering coefficients:")
for node, coeff in top_10_clustering_nodes:
    print(f"Node: {node}, Clustering Coefficient: {coeff}")

# Create subgraph for top 10 nodes by degree centrality
top_degree_nodes_list = [node for node, _ in top_degree_nodes]
subgraph_degree = G.subgraph(top_degree_nodes_list)

# Create subgraph for top 10 nodes by betweenness centrality
top_betweenness_nodes_list = [node for node, _ in top_betweenness_nodes]
subgraph_betweenness = G.subgraph(top_betweenness_nodes_list)

# plot subgrapgs
fig, axes = plt.subplots(1, 2, figsize=(16, 8))
plt.sca(axes[0])
nx.draw_networkx(subgraph_degree, node_size=300, node_color='skyblue', with_labels=True)
plt.title("Top 10 Nodes by Degree Centrality")
plt.sca(axes[1])
nx.draw_networkx(subgraph_betweenness, node_size=300, node_color='lightgreen', with_labels=True)
plt.title("Top 10 Nodes by Betweenness Centrality")
plt.tight_layout()
plt.show()

# Define the interaction map based on top 10 nodes by degree centrality
interaction_map = {
    "Itgax": ["Fcgr3", "Csf1r", "Ptprc", "Cd86"],
    "Fcgr3": ["Itgax", "Tyrobp", "Csf1r", "Ptprc", "Cd86"],
    "Tyrobp": ["Fcgr3", "Csf1r"],
    "Csf1r": ["Tyrobp", "Fcgr3", "Ptprc", "Itgax", "Stat1"],
    "Ptprc": ["Itgax", "Fcgr3", "Csf1r", "Cd86", "Cxcl10"],
    "Cd86": ["Itgax", "Fcgr3", "Csf1r", "Ptprc", "Cxcl10"],
    "Cxcl10": ["Cd86", "Ptprc", "Stat1", "Ifit3", "Irf7"],
    "Stat1": ["Csf1r", "Cxcl10", "Ifit3", "Irf7"],
    "Ifit3": ["Cxcl10", "Stat1", "Irf7"],
    "Irf7": ["Stat1", "Cxcl10", "Ifit3"]}

# Assign interaction strengths (k_ij) arbitrarily for demonstration
interaction_strength = 0.1

# List of proteins
proteins = list(interaction_map.keys())

# Initialize interaction matrix
n_proteins = len(proteins)
interaction_matrix = np.zeros((n_proteins, n_proteins))

# Fill interaction matrix based on interaction map
for i, protein in enumerate(proteins):
    for target in interaction_map[protein]:
        j = proteins.index(target)
        interaction_matrix[i, j] = interaction_strength

# Define the ODE system
def protein_ode(t, concentrations):
    dPdt = np.zeros(n_proteins)
    for i in range(n_proteins):
        # Inputs to the protein i
        input_flux = np.sum(interaction_matrix[:, i] * concentrations)
        # Outputs from the protein i
        output_flux = np.sum(interaction_matrix[i, :] * concentrations[i])
        # Net change in protein concentration
        dPdt[i] = input_flux - output_flux
    return dPdt

# Set initial concentrations: half = 1, one-quarter = 2, one-quarter = 3
initial_concentrations = np.ones(n_proteins)
n_half = n_proteins // 2
n_quarter = n_proteins // 4
initial_concentrations[n_half:n_half + n_quarter] = 2
initial_concentrations[n_half + n_quarter:] = 3

# Time span for the simulation (e.g., 0 to 100 units of time)
t_span = (0, 100)
t_eval = np.linspace(t_span[0], t_span[1], 500)

# Solve the ODE system
solution = solve_ivp(protein_ode, t_span, initial_concentrations, t_eval=t_eval, method='RK45')

# Plot the results
plt.figure(figsize=(12, 8))
for i, protein in enumerate(proteins):
    plt.plot(solution.t, solution.y[i], label=protein) 
plt.title("Protein Concentration Dynamics Over Time")
plt.xlabel("Time")
plt.ylabel("Protein Concentration")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

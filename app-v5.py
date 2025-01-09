import streamlit as st
import pandas as pd
import numpy as np
# group photo
st.image("group-photo.jpeg", caption="Jacobson Lab - UCSF", use_container_width=True)


# Load the primary dataset
@st.cache_data
def load_main_data():
    data_raw = pd.read_csv("HC_clust_output.csv")
    data = data_raw[['gene_name', 'resno_y', 'org_res', 'Count', 
                     'Gain-of-function', 'Inconclusive', 'Likely Gain-of-function', 
                     'Likely Loss-of-function', 'Likely Neutral', 'Loss-of-function', 
                     'Neutral', 'cluster_labels', 'model_no', 'x', 'y', 'z', 'residue_position_id']]
    return data

# Load the variant-level data
@st.cache_data
def load_variant_data():
    variant_data = pd.read_csv("combined-individual-mutations.csv")
    return variant_data[['Legacy Mutation ID', 'gene_name', 'resno', 'org_res', 'mut_res', 'Count']]

# Function to find nearby residues
def find_nearby_residues(data, residue_position_id, distance_cutoff=3.0):
    target_residue = data[data['residue_position_id'] == residue_position_id]
    if target_residue.empty:
        return pd.DataFrame()  # Return empty DataFrame if residue not found

    target_coords = target_residue[['x', 'y', 'z']].values[0]
    data['Distance'] = np.linalg.norm(data[['x', 'y', 'z']].values - target_coords, axis=1)

    nearby = data[(data['Distance'] <= distance_cutoff) & (data['residue_position_id'] != residue_position_id)]
    return nearby[['gene_name', 'resno_y', 'org_res', 'Distance','Count','cluster_labels']].sort_values(by='Distance')

# Load data
main_data = load_main_data()
variant_data = load_variant_data()

# Title
st.title("Jacobson Lab UCSF :Kinase Mutation Explorer")

#st.title("Jacobson lab UCSF :Kinase Atlas")  # Large title
#st.header("Introduction")            # Medium header
#st.subheader("Purpose")              # Small header
#st.text("This is plain text.")       # Plain text
#st.write("This supports **markdown** and _italic_ formatting.")  # Markdown


# Dropdowns for selection
kinases = main_data['gene_name'].unique()
selected_kinase = st.selectbox("Select Kinase Name", kinases, key="mut_kinase")

# Dynamically filter residues based on kinase
filtered_residues = main_data[main_data['gene_name'] == selected_kinase]['resno_y'].unique()
selected_residue = st.selectbox("Select Residue Position", filtered_residues, key="mut_residue")

# Section 1: Mutation Details
st.header("Mutation Details")
filtered_data = main_data[(main_data['gene_name'] == selected_kinase) & (main_data['resno_y'] == selected_residue)]
st.subheader("Mutation Details")
st.write(filtered_data)

# Section 2: Find Nearby Residues
st.header("Find Nearby Residue Positions with Mutations")
nearby_residue_id = main_data[(main_data['gene_name'] == selected_kinase) & (main_data['resno_y'] == selected_residue)]['residue_position_id'].values[0]
nearby_results = find_nearby_residues(main_data, nearby_residue_id)
st.subheader("Nearby Residues")
st.write(nearby_results)

# Download button for nearby residues
st.download_button(
    label="Download Nearby Residues as CSV",
    data=nearby_results.to_csv(index=False).encode('utf-8'),
    file_name="nearby_residues.csv",
    mime="text/csv"
)

# Section 3: Variant-Level Information
st.header("Variant-Level Information")
variant_filtered = variant_data[(variant_data['gene_name'] == selected_kinase) & (variant_data['resno'] == selected_residue)]
st.subheader("Variant-Level Details")
st.write(variant_filtered)

# Download button for variant-level information
st.download_button(
    label="Download Variant Data as CSV",
    data=variant_filtered.to_csv(index=False).encode('utf-8'),
    file_name="variant_level_data.csv",
    mime="text/csv"
)

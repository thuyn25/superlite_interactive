'''
author:        thuyn25 <thuyn27@lsu.edu>
date:          2026-03-07 16:29:36

    This is the main file for the superlite application. It serves as the entry point for the application and is responsible for creating interactive interface that allows users to input their spectra data, deredshift it, and visualize the results with line identification of major relevant species in supernovae spectra.
'''

import streamlit as st
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from astropy import units as u
from specutils import Spectrum1D
import csv

import streamlit as st
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import plotly.express as px

# --- 1. Data Loading --- 
@st.cache_data
def load_line_db():
    try:
        df = pd.read_csv('template/spectral_data.csv')
        # Filter out NaN values for each species column immediately
        spectral_lines = {col: df[col].dropna().tolist() for col in df.columns if col != 'index'}
        major_species = list(spectral_lines.keys())
        return spectral_lines, major_species
    except FileNotFoundError:
        st.error("Line database file not found.")
        return {}, []

LINE_DB, MAJOR_SPECIES = load_line_db()

# Create a stable color mapping for all species
colors = px.colors.qualitative.Alphabet + px.colors.qualitative.Dark24
COLOR_MAP = {species: colors[i % len(colors)] for i, species in enumerate(MAJOR_SPECIES)}

# --- 2. SIDEBAR CONTROLS ---
st.sidebar.title('Controls')
uploaded_file = st.sidebar.file_uploader("Upload Spectrum (lambda, flux)", type=["txt", "csv", "dat"])

st.sidebar.markdown("### Global Shift")
# Added a global redshift to de-redshift the uploaded data
global_z = st.sidebar.number_input("Global Redshift (z)", value=0.0, step=0.001, format="%.4f")

st.sidebar.markdown("### Line Identification")
col1, col2 = st.sidebar.columns(2)
line_params = {}

for i, species in enumerate(MAJOR_SPECIES):
    target_col = col1 if i % 2 == 0 else col2
    with target_col:
        # Display the species name in its assigned color
        st.markdown(f"<span style='color:{COLOR_MAP[species]};font-weight:bold'>{species}</span>", unsafe_allow_html=True)
        is_selected = st.checkbox("Show", key=f"check_{species}", value=(species == "H I" or species == "H"))
        
        if is_selected:
            z_val = st.number_input(f'z ({species})', value=global_z, step=0.01, format="%.4f", key=f'z_{species}', label_visibility="collapsed")
            v_val = st.number_input(f'v_exp ({species})', value=0, step=1000, format="%d", key=f'v_{species}', label_visibility="collapsed")
            line_params[species] = {'z': z_val, 'v_exp': v_val}
        st.markdown("---")

def load_data(file, z):
    data = np.loadtxt(file)
    # Primary axis will be Rest Wavelength
    wavelength_rest = data[:, 0] / (1 + z)
    flux = data[:, 1]
    # Return both rest and original observed wavelengths for axis syncing
    return wavelength_rest, data[:, 0], flux

# --- 4. MAIN INTERFACE ---
st.title("🚀 Supernova Spectroscopic Analyzer")

if uploaded_file:
    wav_rest, wav_obs, flux = load_data(uploaded_file, global_z)

    fig = go.Figure()

    # Plot the main spectrum
    fig.add_trace(go.Scatter(
        x=wav_rest, 
        y=flux, 
        name="Spectrum", 
        line=dict(color='white', width=1)
    ))

    vc = 299792.458 # km/s
    
    # Plot identified lines
    for species, params in line_params.items():
        spec_color = COLOR_MAP[species]
        lines = LINE_DB[species]
        
        # Calculate relative shift: species is shifted by its v_exp and its z relative to the rest frame
        # If species_z == global_z, only v_exp affects the position on the rest axis
        z_factor = (1 + params['z']) / (1 + global_z)
        
        for idx, rest_wav in enumerate(lines):
            shifted_wav = rest_wav * z_factor * (1 - params['v_exp'] / vc)
            
            fig.add_vline(
                x=shifted_wav, 
                line_dash="dot", 
                line_color=spec_color, 
                line_width=1.5,
                opacity=0.8
            )
            
            # Add annotation for the species name
            if idx == 0:
                fig.add_annotation(
                    x=shifted_wav, y=1, yref="paper",
                    text=species, font=dict(color=spec_color),
                    showarrow=False, textangle=-90, yanchor="top"
                )

    # --- DUAL AXIS CONFIGURATION ---
    fig.update_layout(
        template="plotly_dark",
        height=700,
        margin=dict(t=100),
        xaxis=dict(
            title="Rest Wavelength (Å)",
            range=[wav_rest.min(), wav_rest.max()]
        ),
        # Secondary axis (Observed Frame)
        xaxis2=dict(
            title="Observed Wavelength (Å)",
            overlaying='x',
            side='top',
            # This range reflects the actual values in the file
            range=[wav_obs.min(), wav_obs.max()],
            showgrid=False
        ),
        yaxis=dict(title="Flux"),
        hovermode="x unified"
    )

    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("Upload a spectrum to begin.")

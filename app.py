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
import plotly.express as px
st.set_page_config(page_title="🚀 Supernova Spectroscopic Analyzer", layout="wide")

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
uploaded_file = st.sidebar.file_uploader(r"Upload Spectrum ($\lambda$, Flux)", type=["txt", "csv", "dat"])

# st.sidebar.markdown("### Global Shift")
# Added a global redshift to de-redshift the uploaded data
global_z = st.sidebar.number_input("Global Redshift (z)", value=0.0, step=0.001, format="%.4f")

st.sidebar.markdown("### Line Identification")
col1, col2 = st.sidebar.columns(2)
line_params = {}

for i, species in enumerate(MAJOR_SPECIES):
    target_col = col1 if i % 2 == 0 else col2
    with target_col:
        color = COLOR_MAP[species]
        label_html = f"""
            <div style="display: flex; align-items: center; margin-bottom: -10px;">
                <div style="background-color: {color}; width: 6px; height: 18px; margin-right: 8px; border-radius: 2px;"></div>
                <span style="font-weight: bold; font-size: 14px;">{species}</span>
            </div>
        """

        st.markdown(label_html, unsafe_allow_html=True)
        is_selected = st.checkbox("", key=f"check_{species}", value=(species in ["H", "H I"]))
        
        # Display the species name in its assigned color
        # st.markdown(f"<span style='color:{COLOR_MAP[species]};font-weight:bold'>{species}</span>", unsafe_allow_html=True)

        # is_selected = st.checkbox("", key=f"check_{species}", value=(species == "H I" or species == "H"))
        
        if is_selected:
            z_val = st.number_input(f'z = ', value=global_z, step=0.01, format="%.4f", key=f'z_{species}')
            v_val = st.number_input(r'$v_{exp}$ = ', value=0, step=1000, format="%d", key=f'v_{species}')
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

    y_min, y_max = flux.min(), flux.max()
    y_padding = (y_max - y_min) * 0.15

    fig = go.Figure()

    # Plot the main spectrum
    fig.add_trace(go.Scatter(
        x=wav_rest, 
        y=flux, 
        name="Spectrum", 
        showlegend=False,
        line=dict(color='white', width=1)
    ))
    fig.update_yaxes(showticklabels=False)  # Hide y-axis labels for cleaner look
    
    # 2. Ghost Trace (Forces Top Axis to appear)
    fig.add_trace(go.Scatter(
        x=wav_obs, 
        y=flux, 
        xaxis='x2', # THIS IS THE KEY: it links to your xaxis2 definition
        name="Observed Frame",
        marker=dict(opacity=0), # Makes it invisible
        line=dict(width=0),     # No line
        showlegend=False,       # Don't show in legend
        hoverinfo='none'        # Don't show in hover
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
                    x=shifted_wav, y=1+0.00, yref="paper",
                    text=species, font=dict(color=spec_color, size=18),
                    showarrow=True, textangle=-90, yanchor="top", 
                )

    top_min = wav_rest.min() * (1 + global_z)
    top_max = wav_rest.max() * (1 + global_z)

    # --- DUAL AXIS CONFIGURATION ---
    fig.update_layout(
        template="plotly_dark",
        height=700,
        margin=dict(t=100, b=40, l=10, r=10),
        
        xaxis=dict(
            title=dict(text="Rest Wavelength (Å)", font=dict(size=18), standoff=10),
            tickfont=dict(size=14), # Size of the numbers on the axis
            range=[wav_rest.min(), wav_rest.max()],
            automargin=True,
            # minor=dict(
            #     ticklen=4, 
            #     tickcolor="rgba(255,255,255,0.5)", 
            #     showgrid=False
            # ),
            # ticks="inside", # Ensures major ticks are visible
            # tickwidth=1.5
        ),
        
        # --- Secondary X Axis (Top) ---
        xaxis2=dict(
            title=dict(text="Observed Wavelength (Å)", font=dict(size=18), standoff=10),
            tickfont=dict(size=14),
            overlaying='x',
            side='top',
            range=[top_min, top_max],
            showgrid=False,
            # Add Minor Ticks here too
            # minor=dict(
            #     ticklen=10, 
            #     tickcolor="rgba(255,255,255,0.5)", 
            #     showgrid=False
            # ),
            # ticks="inside"
        ),
        yaxis=dict(title=dict(text="Flux", font=dict(size=18)), tickfont=dict(size=14),
                   range=[y_min, y_max + y_padding], automargin=True),
        hovermode="x unified"
    )

    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("Upload a spectrum to begin.")

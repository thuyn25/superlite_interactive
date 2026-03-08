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
st.set_page_config(page_title="SN Spec Analyzer", layout="wide")

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

SPECTRA_COLORS = px.colors.qualitative.Set1

# --- 2. SIDEBAR CONTROLS ---
st.sidebar.title('Controls')
uploaded_files = st.sidebar.file_uploader(r"Upload Spectrum ($\lambda$, Flux)", type=["txt", "csv", "dat"], accept_multiple_files=True)

if len(uploaded_files) > 5:
    st.sidebar.warning('Only the first 5 files will be processed.')
    uploaded_files = uploaded_files[:5]

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
    try:
        data = np.loadtxt(file, usecols=(0, 1), comments='#', unpack=True)
        # Primary axis will be Rest Wavelength
        wavelength_rest = data[0] / (1 + z)
        flux = data[1]

        flux_norm = flux / np.median(flux)  # Normalize flux for better visualization
        # Return both rest and original observed wavelengths for axis syncing
        return wavelength_rest, data[0], flux_norm
    except Exception as e:
        st.error(f"Error loading {file.name}: {e}")
        return None, None, None
# --- 4. MAIN INTERFACE ---
st.title("🚀 Supernova Spectroscopic Analyzer")

if uploaded_files:
    fig = go.Figure()

    all_wav_rest, all_wav_obs= [], []

    for idx, file in enumerate(uploaded_files):
        wav_rest, wav_obs, flux = load_data(file, global_z)

        spec_color = SPECTRA_COLORS[idx % len(SPECTRA_COLORS)]
        # Plot the main spectrum
        fig.add_trace(go.Scatter(
            x=wav_rest, 
            y=flux, 
            name=file.name, 
            showlegend=True,
            line=dict(color=spec_color, width=1.3)
        ))
        fig.update_yaxes(showticklabels=False)  # Hide y-axis labels for cleaner look

        all_wav_rest.extend([wav_rest.min(), wav_rest.max()])
        all_wav_obs.extend([wav_obs.min(), wav_obs.max()])
        # all_flux.extend([flux.min(), flux.max()])
    
        # 2. Ghost Trace (Forces Top Axis to appear)
        fig.add_trace(go.Scatter(
            x=wav_obs, 
            y=flux, 
            xaxis='x2', # THIS IS THE KEY: it links to your xaxis2 definition
            marker=dict(opacity=0), # Makes it invisible
            line=dict(width=0),     # No line
            showlegend=False,       # Don't show in legend
            hoverinfo='none'        # Don't show in hover
        ))

    # y_min, y_max = min(all_flux), max(all_flux)
    # y_padding = (y_max - y_min) * 0.2  # Add
    global_rest_min, global_rest_max = min(all_wav_rest), max(all_wav_rest)
    global_obs_min, global_obs_max = min(all_wav_obs), max(all_wav_obs)

    vc = 299792.458 # km/s

    # Plot identified lines
    for species, params in line_params.items():
        spec_color = COLOR_MAP[species]
        lines = LINE_DB[species]
        
        # Calculate relative shift: species is shifted by its v_exp and its z relative to the rest frame
        # If species_z == global_z, only v_exp affects the position on the rest axis
        z_factor = (1 + params['z']) / (1 + global_z)
        
        for l_idx, rest_wav in enumerate(lines):
            shifted_wav = rest_wav * z_factor * (1 - params['v_exp'] / vc)
            
            fig.add_vline(
                x=shifted_wav, 
                line_dash="dot", 
                line_color=spec_color, 
                line_width=1.5,
                opacity=0.8
            )
            
            # Add annotation for the species name
            if l_idx == 0:
                fig.add_annotation(
                    x=shifted_wav, y=1+0.00, yref="paper",
                    text=species, font=dict(color=spec_color, size=18),
                    showarrow=True, textangle=-90, yanchor="top", 
                )

    # --- DUAL AXIS CONFIGURATION ---
    fig.update_layout(
        template="plotly_dark",
        height=750,
        margin=dict(t=100, b=40, l=10, r=10),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, font=dict(size=14)),
        
        xaxis=dict(
            title=dict(text="Rest Wavelength (Å)", font=dict(size=18), standoff=10),
            tickfont=dict(size=14), # Size of the numbers on the axis
            range=[2000, 10000],
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
            # range=[global_obs_min, global_obs_max],
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
                automargin=True),
        hovermode="x unified"
    )

    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("Upload a spectrum to begin.")

'''
author:        Trang Huynh <thuyn27@lsu.edu>
date:          2026-03-07 15:14:10

    Spectra references come from 3 sources:
        'spectra_harvard.txt': [https://lweb.cfa.harvard.edu/amp/ampdata/kurucz23/sekur.html]
        'spectra_atomicnet.csv': [https://www.atomic-spectra.net]
        'spectra_NIST.csv': [https://physics.nist.gov/PhysRefData/ASD/lines_form.html]

'''

import re
import numpy as np
import pandas as pd
import roman

def strip_fotmatting(val):  # for NIST database
    if isinstance(val, str):
        # This regex removes the leading =" and the trailing "
        # It handles cases like ="2000.0225" -> 2000.0225
        # and ="" -> empty string
        val = re.sub(r'^="', '', val)
        val = re.sub(r'"$', '', val)
    return val

# Loading database #
wavelengths, isotopes = [], []
with open('database/spectra_harvard.txt', 'r') as f:
    for line in f:
        parts = line.split()
        try:
            wl = float(parts[0])

            el_name = f"{parts[3]} {parts[4]}"

            wavelengths.append(wl)
            isotopes.append(el_name)
        except (ValueError, IndexError): continue

wavelengths = np.array(wavelengths)
isotopes = np.array(isotopes)

with open('database/spectra_harvard.csv', 'w') as f:
    f.write('Element,Wl_nm\n')
    for elem, wav in zip(isotopes, wavelengths):
        f.write(f'{elem},{wav}\n')

df_harvard = pd.read_csv('database/spectra_harvard.csv', delimiter=',')
df_atomicnet = pd.read_csv('database/spectra_atomicnet.csv')
NIST_df2 = pd.read_csv('database/spectra_NIST.csv', delimiter=',')
df_NIST =  NIST_df2.applymap(strip_fotmatting)

'''
    Elements to extract for supernovae' spectral lines: 
    ['H I', 'He I', 'He II', 'C I','C II', 'C III', 'C IV', 'N I','N II', 'N III', 'N IV', 'N V', 'O I','[O I]', 'O II','[O II]','[O III]', 'O V', 'O VI','Na I', 'Mg I', 'Mg II', 'Si II', 'S I', 'S II', 'Ca II', '[Ca II]', 'Fe II', 'Fe I']
'''
f1 = open(f'template/spectra_template.csv', 'w')

needed_isotopes = ['H I', 'He I', 'He II', 'C I','C II', 'C III', 'C IV', 'N I','N II', 'N III', 'N IV', 'N V', 'O I', 'O II', 'O V', 'O VI','Na I', 'Mg I', 'Mg II', 'Si II', 'S I', 'S II', 'Ca II', 'Fe II', 'Fe I']
forbidden_isotopes = ['O I', 'O II', 'O III', 'Ca II']
# Extracting from harvard file #
f1.write('Element,Wavelength_A\n')
remaining_isotopes = []
for iso in needed_isotopes:
    if iso in df_harvard['Element'].unique():
        # print(f'iso: {iso} in the list.')
        wav = df_harvard[df_harvard['Element'] == iso]['Wl_nm'].values
        for w in wav:
            f1.write(f'{iso},{w * 10}\n')
            # f1.write(f'{iso},{w}\n')
    elif iso not in df_harvard['Element'].unique():
        # print(f'iso: {iso} not in the list.')
        remaining_isotopes.append(iso)

# Extracting from atomicnet file #
# print(f'Isotopes not existing in harvard file: {remaining_isotopes}\n')
# f1.write('atomicnet separating line...\n')
remaining_isotopes2 = []
for iso in remaining_isotopes:
    parts = iso.split()
    el = parts[0]
    ion = parts[1]
    df = df_atomicnet[(df_atomicnet['Element'] == el) & (df_atomicnet['Spectrum'] == ion)]

    if not df.empty:
        wav = df['Wavelength, A'].values
        for w in wav:
            if w >= 2000.0 and w <= 10000.0:
                f1.write(f'{iso},{w}\n')
    elif df.empty:
        remaining_isotopes2.append(iso)

# print(f'Isotopes not existing in atomicnet file: {remaining_isotopes2}\n')

# Extracting from NIST file for forbidden isotopes #
# print(remaining_isotopes)
remaining_forbidden_isotopes = []
print('Writing forbidden lines...')
# f1.write('NIST separating line...\n')
forbidden_df = df_NIST[df_NIST['Type'].isin(['E2', 'M1'])]

for iso in forbidden_isotopes:
    parts = iso.split()
    el = parts[0]
    ion = parts[1]
    ion = roman.fromRoman(ion)
    # print(ion, type(ion))
    df = forbidden_df[(forbidden_df['element'] == el) & (forbidden_df['sp_num'] == ion)]
    if not df.empty:
        wav = df['obs_wl_air(A)'].values
        for w in wav:
            w = float(w)
            if w >= 2000.0 and w <= 10000.0:
                f1.write(f'[{iso}],{w}\n')
    elif df.empty:
        remaining_forbidden_isotopes.append(iso)

print(f'Isotopes not found in any database: {remaining_isotopes2}, [{remaining_forbidden_isotopes}]\n')

f1.close()
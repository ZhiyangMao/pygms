"""
Build the NGA-Sub ground motion databases (interface and intraslab) for pygms.

Inputs (public NGA-Sub products):
- NGAsub_MegaFlatfile_RotD50_050_R211022_public.xlsx
- NGAsubDb_PSA_H1_050_public_R230104.csv, NGAsubDb_PSA_H2_050_public_R230104.csv
- NGAsub_DurationData_R230609_Public.csv

Records are classified by the flatfile column Intra_Inter_Flag (0: interface, 1: intraslab).
H1 and H2 components are stored as separate records (H1 block followed by H2 block).
Sa is interpolated to the pygms period grid (0.01 s to 10 s, 0.01 s step); missing values
(negative in the NGA-Sub tables) are set to 1e-15. Sa is written with 5 significant digits
into xz-compressed chunks.

Usage:
    python build_gmdb_nga_sub.py --raw-dir <folder with the NGA-Sub files> --out-dir ../gmdb
"""
import os
import io
import lzma
import argparse
import numpy as np
import pandas as pd

FLATFILE = 'NGAsub_MegaFlatfile_RotD50_050_R211022_public.xlsx'
PSA_H1 = 'NGAsubDb_PSA_tables_050_public_R230104/NGAsubDb_PSA_H1_050_public_R230104.csv'
PSA_H2 = 'NGAsubDb_PSA_tables_050_public_R230104/NGAsubDb_PSA_H2_050_public_R230104.csv'
DURATION = 'NGAsub_DurationData_R230609_Public.csv'

EVENT_TYPES = {'interface': 0, 'intraslab': 1}
PERIODS = np.round(np.arange(0.01, 10.01, 0.01), 2)
CHUNK_ROWS = 20000


def parse_period(col):
    # e.g., 'T0pt010S' -> 0.01, 'T10pt000S' -> 10.0
    return float(col[1:-1].replace('pt', '.'))


def align(df, rsn, name):
    # reorder df to follow rsn exactly and fail loudly if any record is missing
    out = df.set_index('NGAsubRSN').reindex(rsn)
    if out.index.has_duplicates or out.isnull().all(axis=1).any():
        raise ValueError('{}: records missing or duplicated for the selected RSNs.'.format(name))
    return out.reset_index()


def interp_psa(psa):
    cols = [c for c in psa.columns if c.startswith('T') and c.endswith('S') and 'pt' in c]
    t = np.array([parse_period(c) for c in cols])
    order = np.argsort(t)
    t, vals = t[order], psa[cols].to_numpy(dtype=float)[:, order]
    out = np.array([np.interp(PERIODS, t, row) for row in vals])
    out[out < 0] = 1e-15
    return out


def write_psa_chunks(psa, out_dir):
    n_chunk = int(np.ceil(len(psa) / CHUNK_ROWS))
    for i in range(n_chunk):
        buf = io.StringIO()
        pd.DataFrame(psa[i*CHUNK_ROWS:(i+1)*CHUNK_ROWS]).to_csv(buf, header=False, index=False, float_format='%.5g')
        with lzma.open(os.path.join(out_dir, 'psa_p{}.csv.xz'.format(i+1)), 'wb', preset=6) as f:
            f.write(buf.getvalue().encode())
    return n_chunk


def build(raw_dir, out_dir, event_type, flatfile, psa_h1, psa_h2, duration):
    ff = flatfile[flatfile['Intra_Inter_Flag'] == EVENT_TYPES[event_type]].reset_index(drop=True)
    rsn = ff['NGAsubRSN']
    h1 = align(psa_h1, rsn, 'PSA H1')
    h2 = align(psa_h2, rsn, 'PSA H2')
    du = align(duration, rsn, 'Duration')

    cur_dir = os.path.join(out_dir, 'nga_sub_{}'.format(event_type))
    os.makedirs(cur_dir, exist_ok=True)

    # file names: NGAsubRSN<rsn>_<component file name>
    names = []
    for df, col in [(h1, 'accFilePathH1'), (h2, 'accFilePathH2')]:
        base = df[col].astype(str).str.replace('\\', '/', regex=False).str.split('/').str[-1]
        names.append('NGAsubRSN' + rsn.astype(str) + '_' + base)
    tables = {
        'filename.csv': pd.concat(names, ignore_index=True),
        'Earthquake Magnitude.csv': pd.concat([ff['Earthquake_Magnitude']]*2, ignore_index=True),
        'EpiD (km).csv': pd.concat([ff['EpiD_km']]*2, ignore_index=True),
        'Preferred Vs30 (m-s).csv': pd.concat([ff['Vs30_Selected_for_Analysis_m_s']]*2, ignore_index=True),
        'ds575.csv': pd.concat([du['H1_Duration_0pt05_0pt75'], du['H2_Duration_0pt05_0pt75']], ignore_index=True),
        'ds595.csv': pd.concat([du['H1_Duration_0pt05_0pt95'], du['H2_Duration_0pt05_0pt95']], ignore_index=True),
    }
    for fname, ser in tables.items():
        ser.to_csv(os.path.join(cur_dir, fname), index=False, header=False)
    pd.DataFrame([PERIODS]).to_csv(os.path.join(cur_dir, 'period.csv'), index=False, header=False)

    psa = np.vstack([interp_psa(h1), interp_psa(h2)])
    n_chunk = write_psa_chunks(psa, cur_dir)
    print('{}: {} records ({} rows with H1 and H2), {} Sa chunks -> {}'.format(
        event_type, len(rsn), len(psa), n_chunk, cur_dir))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-dir', required=True)
    parser.add_argument('--out-dir', default=os.path.join(os.path.dirname(__file__), '..', 'gmdb'))
    parser.add_argument('--event-types', nargs='+', default=list(EVENT_TYPES.keys()))
    args = parser.parse_args()

    flatfile = pd.read_excel(os.path.join(args.raw_dir, FLATFILE))
    psa_h1 = pd.read_csv(os.path.join(args.raw_dir, PSA_H1))
    psa_h2 = pd.read_csv(os.path.join(args.raw_dir, PSA_H2))
    duration = pd.read_csv(os.path.join(args.raw_dir, DURATION))
    for cur_type in args.event_types:
        build(args.raw_dir, args.out_dir, cur_type, flatfile, psa_h1, psa_h2, duration)

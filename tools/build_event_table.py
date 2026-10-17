"""
Build the earthquake event tables used by the EventPenalty option in pygms.

Each table is row-aligned with the database it belongs to and has the columns
Filename, RSN, EQID, EarthquakeName, taken from the NGA-West2 flatfile.

Outputs:
- gmdb/nga_west_event.csv (default database; aligned with gmdb/nga_west_filename.csv)
- gmdb/nga_west2/event.csv (NGA-West2 database; aligned with gmdb/nga_west2/rsn.csv)

Usage:
    python build_event_table.py --flatfile Updated_NGA_West2_Flatfile_RotD50_d050_public_version.xlsx
"""
import os
import argparse
import pandas as pd

GMDB = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gmdb')


def build(flatfile, filenames, rsn):
    ff = flatfile[['Record Sequence Number', 'EQID', 'Earthquake Name']].set_index('Record Sequence Number')
    if not rsn.isin(ff.index).all():
        raise ValueError('Some RSNs are not found in the flatfile.')
    ff = ff.loc[rsn.values]
    return pd.DataFrame({
        'Filename': filenames.values,
        'RSN': rsn.values,
        'EQID': ff['EQID'].astype(int).map('{:04d}'.format).values,
        'EarthquakeName': ff['Earthquake Name'].astype(str).str.strip().values,
    })


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--flatfile', required=True)
    args = parser.parse_args()
    flatfile = pd.read_excel(args.flatfile)

    # default database: file names such as RSN1_HELENA.A_A-HMC180.txt
    filenames = pd.read_csv(os.path.join(GMDB, 'nga_west_filename.csv'), header=None)[0]
    rsn = filenames.str.extract(r'^RSN(\d+)_')[0].astype(int)
    build(flatfile, filenames, rsn).to_csv(os.path.join(GMDB, 'nga_west_event.csv'), index=False)

    # NGA-West2 database: one record per RSN
    rsn = pd.read_csv(os.path.join(GMDB, 'nga_west2', 'rsn.csv'), header=None)[0].astype(int)
    build(flatfile, rsn, rsn).to_csv(os.path.join(GMDB, 'nga_west2', 'event.csv'), index=False)

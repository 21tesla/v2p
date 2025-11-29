import sqlite3
import genomicsqlite
import sys
import pandas as pd
import os

PROJECT_DIR = os.environ['V2P_DIR']
DATABASE_LOCATION = sys.argv[3]
NUCS = ['A','T','G','C', 'a', 't', 'g', 'c']

CUTOFFS = {'HP:0033127': 0.1567193267520465,
 'HP:0040064': 0.1462383623540305,
 'HP:0000707': 0.218444388264641,
 'HP:0001939': 0.1343807652765258,
 'HP:0000152': 0.2044814790155267,
 'HP:0001626': 0.1778303343086918,
 'HP:0000119': 0.1660729258282984,
 'HP:0000478': 0.3116954769244471,
 'HP:0002715': 0.1708217552071461,
 'HP:0001574': 0.1733027692065684,
 'HP:0001871': 0.1351526870849739,
 'HP:0025031': 0.1604417003743339,
 'HP:0002664': 0.2689097200492832,
 'HP:0002086': 0.0989268196827497,
 'HP:0000818': 0.0525938291132373,
 'HP:0000598': 0.2912918120493089,
 'HP:0025354': 0.0876004261284816,
 'HP:0001197': 0.0030211178500946,
 'HP:0001507': 0.0330074809050576,
 'HP:0025142': 0.0398178148872888,
 'HP:0000769': 0.133028045951846,
 'HP:0001608': 0.003593885070548,
 'HP:0045027': 0.1,
 'Pathogenic': 0.3985216152439556}

 
NAMES_HPO = {'Musculoskeletal': 'HP:0033127',
 'Limbs': 'HP:0040064',
 'Nervous': 'HP:0000707',
 'Metabolism/homeostasis': 'HP:0001939',
 'Head/neck': 'HP:0000152',
 'Cardiovascular': 'HP:0001626',
 'Genitourinary': 'HP:0000119',
 'Eye': 'HP:0000478',
 'Immune': 'HP:0002715',
 'Integument': 'HP:0001574',
 'Blood/blood-forming tissues': 'HP:0001871',
 'Digestive': 'HP:0025031',
 'Neoplasm': 'HP:0002664',
 'Respiratory': 'HP:0002086',
 'Endocrine': 'HP:0000818',
 'Ear': 'HP:0000598',
 'Cellular': 'HP:0025354',
 'Prenatal development/birth': 'HP:0001197',
 'Growth': 'HP:0001507',
 'Constitutional': 'HP:0025142',
 'Breast': 'HP:0000769',
 'Voice': 'HP:0001608',
 'Thoracic cavity': 'HP:0045027',
 'Pathogenic': 'Pathogenic'}

ORDER = ['ID', 'Musculoskeletal', 'Limbs', 'Nervous', 'Metabolism/homeostasis', 'Head/neck', 'Cardiovascular', 'Genitourinary', 'Eye',
    'Immune', 'Integument', 'Blood/blood-forming tissues', 'Digestive', 'Neoplasm', 'Respiratory', 'Endocrine', 'Ear', 'Cellular',
    'Prenatal development/birth', 'Growth', 'Constitutional', 'Breast', 'Voice','Thoracic cavity', 'Pathogenic']
OUT_ORDER = ['ID','V2P_predicted_phenotypes', 'Musculoskeletal', 'Limbs', 'Nervous', 'Metabolism/homeostasis', 'Head/neck', 'Cardiovascular', 'Genitourinary', 'Eye',
    'Immune', 'Integument', 'Blood/blood-forming tissues', 'Digestive', 'Neoplasm', 'Respiratory', 'Endocrine', 'Ear', 'Cellular',
    'Prenatal development/birth', 'Growth', 'Constitutional', 'Breast', 'Voice','Thoracic cavity', 'Pathogenic']

data = pd.read_csv(sys.argv[1], sep='\t', low_memory=False)
data['#CHROM'] = data['#CHROM'].apply(lambda v: str(v).replace('chr',''))
data['ID'] = data.apply(lambda r: '_'.join([str(r['#CHROM']), str(int(r['POS'])),r['REF'], r['ALT']]), axis=1) 
chroms = list(data['#CHROM'].unique())

tmpfile = sys.argv[2]

def get_crisp(df):
    crisp = []
    for _, r in df.iterrows():
        rc = [l for p, l in zip(r, df.columns) if l != 'ID' and p >= CUTOFFS[NAMES_HPO[l]]]
        if not rc:
            rc = ['Benign']
        crisp.append(','.join(rc))
    return crisp

def get_predictions(conn, variants):
    cursor = conn.cursor()
    escaped_columns = ','.join([f'"{col}"' for col in ORDER])
    query = f'''
    SELECT {escaped_columns}
    FROM predictions
    WHERE ID IN {repr(tuple(map(str, variants['ID'].tolist() + [''])))}
    '''
    cursor.execute(query)
    return cursor.fetchall()

# Selects precomputed predictions from the downloaded database based on the variant
# Chromosome_Position_Ref_Alt ID
for chrom in chroms:

    chrom_df = data.loc[data['#CHROM'] == chrom]
    snps = chrom_df.loc[((chrom_df['REF'].isin(NUCS)) & (chrom_df['ALT'].isin(NUCS)))]
    indels = chrom_df.loc[~chrom_df['ID'].isin(snps['ID'])]
    
    snp_dbconn = genomicsqlite.connect(
        DATABASE_LOCATION + '/chrom' + str(chrom) + '_compressed.db',
        read_only=True,
    )

    indel_dbconn = genomicsqlite.connect(
        DATABASE_LOCATION + '/gnomad_chrom' + str(chrom) + '_compressed.db',
        read_only=True,
    )

    snp_rows = pd.DataFrame(columns=ORDER)
    indel_rows = pd.DataFrame(columns=ORDER)
    if not snps.empty:
        snp_rows = pd.DataFrame(get_predictions(snp_dbconn, snps), columns=ORDER)
    if not indels.empty:
        indel_rows = pd.DataFrame(get_predictions(indel_dbconn, indels), columns=ORDER)

    snp_rows['V2P_predicted_phenotypes'] = get_crisp(snp_rows)
    indel_rows['V2P_predicted_phenotypes'] = get_crisp(indel_rows)

    file_populated = os.path.exists(tmpfile) and os.stat(tmpfile).st_size > 0
    header = None if file_populated else True
    mode = 'a+' if file_populated else 'w+'
    snp_rows[OUT_ORDER].to_csv(tmpfile, index=None, header=header, mode=mode)
    indel_rows[OUT_ORDER].to_csv(tmpfile, index=None, header=None, mode='a+')

found = pd.read_csv(tmpfile, usecols=['ID'])
data = data.loc[~data['ID'].isin(found['ID'])]
outname = sys.argv[1].split('/')[-1].replace('.vcf', '_novel.vcf')
data.to_csv(outname, index=None, sep='\t')

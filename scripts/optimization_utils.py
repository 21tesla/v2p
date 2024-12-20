import time
import sys
import pandas as pd
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import MinMaxScaler, OrdinalEncoder
import random
from tqdm import tqdm
try:
    from imblearn.pipeline import Pipeline
except:
    from sklearn.pipeline import Pipeline
    print('Warning: imbalanced-learn not available', file=sys.stderr)

NEGONE_FEATURES = ['BLOSUM62','ProteinLengthChange','MaxEntScan_alt','MaxEntScan_diff','MaxEntScan_ref','ada_score','rf_score','rel_cDNA_pos','rel_CDS_pos','rel_prot_pos','Phosphorylation','Acetylation','Methylation','Ubiquitination','Glycosylation','PTM','RSA','ASA','AF_Relative_ASA','IUPRED2','ANCHOR2','A3D_SCORE','n_contacts','distance_com','concavity_score','S_DDG[SEQ]','S_DDG[3D]','hgmd_mutcount','gnomsingle_mutcount','gnom_mutcount','AF_confidence','isHomomultimer','num_interactions','ppi_combined_0','ppi_combined_1','ppi_combined_2','ppi_combined_3','ppi_combined_4','ppi_combined_5','ppi_combined_6','ppi_combined_7','ppi_combined_8','ppi_combined_9','ppi_combined_10','ppi_combined_11','ppi_combined_12','ppi_combined_13','ppi_combined_14','ppi_combined_15','ppi_combined_16','ppi_combined_17','ppi_combined_18','ppi_combined_19','ppi_combined_20','ppi_combined_21','ppi_combined_22','ppi_combined_23','ppi_combined_24','ppi_combined_25','ppi_combined_26','ppi_combined_27','ppi_combined_28','ppi_combined_29','ppi_combined_30','ppi_combined_31','ppi_combined_32','ppi_combined_33','ppi_combined_34','ppi_combined_35','ppi_combined_36','ppi_combined_37','ppi_combined_38','ppi_combined_39','ppi_combined_40','ppi_combined_41','ppi_combined_42','ppi_combined_43','ppi_combined_44','ppi_combined_45','ppi_combined_46','ppi_combined_47','ppi_combined_48','ppi_combined_49','ppi_combined_50','ppi_combined_51','ppi_combined_52','ppi_combined_53','ppi_combined_54','ppi_combined_55','ppi_combined_56','ppi_combined_57','ppi_combined_58','ppi_combined_59','ppi_combined_60','ppi_combined_61','ppi_combined_62','ppi_combined_63','DRNApredDNAscore_aa','ASAquick_normscore_aa','ASAquick_rawscore_aa','DFLpredScore_aa','DRNApredRNAscore_aa','DisoDNAscore_aa','DisoPROscore_aa','DisoRNAscore_aa','MMseq2_conservation_level_aa','MMseq2_conservation_score_aa','MoRFchibiScore_aa','PSIPRED_helix_aa','PSIPRED_strand_aa','SCRIBERscore_aa','SignalP_score_aa','PHOSPHORYLATION','ACETYLATION','UBIQUITINATION','S-NITROSYLATION','N-GLYCOSYLATION','METHYLATION','O-GLYCOSYLATION','MYRISTOYLATION','C-GLYCOSYLATION','SUMOYLATION','S-GLYCOSYLATION','polyphen_nobs','polyphen_normasa','polyphen_dvol','polyphen_dprop','polyphen_bfact','polyphen_hbonds','polyphen_avenhet','polyphen_mindhet','polyphen_avenint','polyphen_mindint','polyphen_avensit','polyphen_mindsit','polyphen_idpmax','polyphen_idpsnp','polyphen_idqmin','MOD_RES','REGION','INTERACTION_REGION','REQUIRED_FOR_INTER','Dst2Splice','Grantham','SpliceAI-acc-gain','SpliceAI-acc-loss','SpliceAI-don-gain','SpliceAI-don-loss','MMSp_acceptorIntron','MMSp_acceptor','MMSp_exon','MMSp_donor','MMSp_donorIntron','dbscSNV-ada_score','dbscSNV-rf_score','Charge','Volume','Hydrophobicity','Polarity','Ex','PAM250','BLOSUM','JM','HGMD2003','VB','Transition','COSMIC','COSMICvsSWISSPROT','HAPMAP','COSMICvsHAPMAP','ATP_binding_gbind','Ca2+_binding_gbind','DNA_binding_gbind','HEME_binding_gbind','Mg2+_binding_gbind','Mn2+_binding_gbind','RNA_binding_gbind']
MEDIAN_FEATURES = ['Conservation','TSSDistance','Eigen_PC_raw_coding','Eigen_raw_coding','GERPplus_plus_NR','GERPplus_plus_RS','GM12878_confidence_value','GM12878_fitCons_score','GenoCanyon_score','H1_hESC_confidence_value','H1_hESC_fitCons_score','HUVEC_confidence_value','HUVEC_fitCons_score','LINSIGHT','LIST_S2_score','LRT_Omega','LRT_score','MPC_score','MutationAssessor_score','SiPhy_29way_logOdds','bStatistic_x','integrated_confidence_value','integrated_fitCons_score','phastCons100way_vertebrate','phastCons17way_primate','phastCons30way_mammalian','phyloP100way_vertebrate','phyloP17way_primate','phyloP30way_mammalian','GDI','MSC_95CI','Selective_pressure','Clarks_distance','CDS_len','Number_of_paralogs','denovo_Zscore','RVIS','Indispensability_score','NearestExonJB_distance','s_het','gtex_Adipose_-_Subcutaneous','gtex_Adipose_-_Visceral_(Omentum)','gtex_Adrenal_Gland','gtex_Artery_-_Aorta','gtex_Artery_-_Coronary','gtex_Artery_-_Tibial','gtex_Bladder','gtex_Brain_-_Amygdala','gtex_Brain_-_Anterior_cingulate_cortex_(BA24)','gtex_Brain_-_Caudate_(basal_ganglia)','gtex_Brain_-_Cerebellar_Hemisphere','gtex_Brain_-_Cerebellum','gtex_Brain_-_Cortex','gtex_Brain_-_Frontal_Cortex_(BA9)','gtex_Brain_-_Hippocampus','gtex_Brain_-_Hypothalamus','gtex_Brain_-_Nucleus_accumbens_(basal_ganglia)','gtex_Brain_-_Putamen_(basal_ganglia)','gtex_Brain_-_Spinal_cord_(cervical_c-1)','gtex_Brain_-_Substantia_nigra','gtex_Breast_-_Mammary_Tissue','gtex_Cells_-_Cultured_fibroblasts','gtex_Cells_-_EBV-transformed_lymphocytes','gtex_Cervix_-_Ectocervix','gtex_Cervix_-_Endocervix','gtex_Colon_-_Sigmoid','gtex_Colon_-_Transverse','gtex_Esophagus_-_Gastroesophageal_Junction','gtex_Esophagus_-_Mucosa','gtex_Esophagus_-_Muscularis','gtex_Fallopian_Tube','gtex_Heart_-_Atrial_Appendage','gtex_Heart_-_Left_Ventricle','gtex_Kidney_-_Cortex','gtex_Kidney_-_Medulla','gtex_Liver','gtex_Lung','gtex_Minor_Salivary_Gland','gtex_Muscle_-_Skeletal','gtex_Nerve_-_Tibial','gtex_Ovary','gtex_Pancreas','gtex_Pituitary','gtex_Prostate','gtex_Skin_-_Not_Sun_Exposed_(Suprapubic)','gtex_Skin_-_Sun_Exposed_(Lower_leg)','gtex_Small_Intestine_-_Terminal_Ileum','gtex_Spleen','gtex_Stomach','gtex_Testis','gtex_Thyroid','gtex_Uterus','gtex_Vagina','gtex_Whole_Blood','haplo','haplo_imputed','GC','CpG','motifECount','motifEHIPos','motifEScoreChng','minDistTSS','minDistTSE','priPhCons','mamPhCons','verPhCons','priPhyloP','mamPhyloP','verPhyloP','targetScan','mirSVR-Score','mirSVR-E','mirSVR-Aln','cHmm_E1','cHmm_E2','cHmm_E3','cHmm_E4','cHmm_E5','cHmm_E6','cHmm_E7','cHmm_E8','cHmm_E9','cHmm_E10','cHmm_E11','cHmm_E12','cHmm_E13','cHmm_E14','cHmm_E15','cHmm_E16','cHmm_E17','cHmm_E18','cHmm_E19','cHmm_E20','cHmm_E21','cHmm_E22','cHmm_E23','cHmm_E24','cHmm_E25','GerpRS','GerpRSpval','GerpN','GerpS','tOverlapMotifs','motifDist','EncodeH3K4me1-sum','EncodeH3K4me1-max','EncodeH3K4me2-sum','EncodeH3K4me2-max','EncodeH3K4me3-sum','EncodeH3K4me3-max','EncodeH3K9ac-sum','EncodeH3K9ac-max','EncodeH3K9me3-sum','EncodeH3K9me3-max','EncodeH3K27ac-sum','EncodeH3K27ac-max','EncodeH3K27me3-sum','EncodeH3K27me3-max','EncodeH3K36me3-sum','EncodeH3K36me3-max','EncodeH3K79me2-sum','EncodeH3K79me2-max','EncodeH4K20me1-sum','EncodeH4K20me1-max','EncodeH2AFZ-sum','EncodeH2AFZ-max','EncodeDNase-sum','EncodeDNase-max','EncodetotalRNA-sum','EncodetotalRNA-max','Dist2Mutation','Freq100bp','Rare100bp','Sngl100bp','Freq1000bp','Rare1000bp','Sngl1000bp','Freq10000bp','Rare10000bp','Sngl10000bp','RemapOverlapTF','RemapOverlapCL']

class DropCorrelatedTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, correlation_cutoff=.95):
        self.correlation_cutoff = correlation_cutoff
        self.to_drop = None

    def fit(self, X, y=None):
        X = pd.DataFrame(X)
        cor_matrix = X.corr().abs()
        upper_tri = cor_matrix.where(np.triu(np.ones(cor_matrix.shape),k=1).astype(np.bool))
        self.to_drop = [column for column in upper_tri.columns if any(upper_tri[column] >= self.correlation_cutoff)]
        return self

    def transform(self, X, y=None):
        X = pd.DataFrame(X)
        X = X.drop(columns=self.to_drop)
        return X.to_numpy()

class SubsetTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, max_columns):
        self.max_columns = max_columns
        self.keep_cols = None

    def fit(self, X, y=None):
        print(X.shape)
        cols = random.sample(list(range(X.shape[1])), round(X.shape[1] * self.max_columns))
        self.keep_cols = cols
        return self

    def transform(self, X, y=None):
        X = pd.DataFrame(X)
        X = X[self.keep_cols]
        print(X.shape)
        return X.to_numpy()

class FillNA(BaseEstimator, TransformerMixin):
    def __init__(self, fill=0):
        self.fill = 0

    def fit(self, X, y=None):
        self.X = X
        return self

    def transform(self, X, y=None):
        self.X = X
        return X.fillna(0)
    
    def get_feature_names_out(self, out):
        return self.X.columns.tolist()

class FastImputer(BaseEstimator, TransformerMixin):
    def __init__(self):
        pass

    def fit(self, X, y=None):
        print('fitting')
        self.medians = np.nanmean(X, axis=0)
        print('done fitting')
        return self

    def transform(self, X, y=None):
        X = X.copy()
        for i in tqdm(range(X.shape[1]), total = X.shape[1]):
            X[:,i][X[:,i] == np.nan] = self.medians[i]
        return X

def drop_allnan(data):
    for col in data.columns:
        if data[col].isna().sum() == len(data):
            data = data.drop(columns=col)
    return data

class BooltoInt(BaseEstimator, TransformerMixin):
    def __init__(self):
        pass

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return X.replace([False, True], [0, 1]) 

def generate_preprocessor(numeric_features, categorical_features, bool_features, N_JOBS, do_categorical=True):

    categorical_transformer = OrdinalEncoder(handle_unknown='use_encoded_value', unknown_value=-1, encoded_missing_value=-1)

    numeric_transformer = Pipeline(steps=[
        ('imputer', FillNA()),
        ('scaler', MinMaxScaler(feature_range =(0, 1), clip=True))])

    transformers = [('numeric', numeric_transformer, numeric_features),
                    ('bool', BooltoInt(), bool_features)]
    if do_categorical:
        transformers.append(('cat', categorical_transformer, categorical_features))
    else:
        transformers.append(('remainder', 'passthrough'))

    preprocessor = ColumnTransformer(
        transformers=[
            ('numeric', numeric_transformer, numeric_features),
            ('cat', categorical_transformer, categorical_features),
        ])

    steps = [('initial', preprocessor)]
    if do_categorical:
        steps.append(('variance_threshold', VarianceThreshold(threshold=0)))
    preprocessor = Pipeline(steps=steps)
    return preprocessor

def preprocess(preprocessor, train_data, train_labels, test_data, quiet=False):
    for k, v in preprocessor.steps:
        if k == 'initial':
            start = time.time()
            v.fit(train_data)
            train_data = pd.DataFrame(v.transform(train_data), columns=v.get_feature_names_out())
            end = time.time()
            if not quiet:
                print(k + ' took ' + str(end - start) + ' to run.')
        elif k == 'oversampling' or k == 'undersampling':
            start = time.time()
            train_data, train_labels = v.fit_resample(train_data, train_labels)
            end = time.time()
            if not quiet:
                print(k + ' took ' + str(end - start) + ' to run.')
        else:
            start = time.time()
            train_data = v.fit_transform(train_data, train_labels)
            end = time.time()
            if not quiet:
                print(k + ' took ' + str(end - start) + ' to run')
    for k, v in preprocessor.steps:
        if k == 'initial':
            test_data = pd.DataFrame(v.transform(test_data), columns=v.get_feature_names_out())
        elif k == 'oversampling' or k == 'undersampling':
            continue
        else:
            test_data = v.transform(test_data)

    train_data = pd.DataFrame(train_data)
    test_data = pd.DataFrame(test_data)

    for col in train_data.columns:
        try:
            train_data[col] = train_data[col].astype('float')
            test_data[col] = test_data[col].astype('float')
        except:
            train_data[col] = train_data[col].astype('category')
            test_data[col] = test_data[col].astype('category')

    return train_data.to_numpy(), train_labels, test_data.to_numpy()

class TimeoutException(Exception):   # Custom exception class
    pass


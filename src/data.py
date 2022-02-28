import os
from glob import glob
from scipy.io.matlab import loadmat
import numpy as np
import pandas as pd
import re

data_root = os.path.abspath(
        os.path.join(
            os.path.dirname(
                os.path.dirname(__file__)
            ),
            'data'
        )
)

class Julich:
    
    ds_root     = os.path.join(data_root, "external", "Julich")
    ds_external = os.path.join(data_root, "external")
    
    def __init__(self):
        self.data_root = data_root
    
    def parcellation_400(self):
        d = os.path.join(self.ds_external)
        separator   = ''
        parce_frame = pd.read_csv(separator.join([d,'/Schaefer2018_400Parcels_7Networks_tab.txt']), delimiter=",")
        parce_list  = parce_frame["Var1"].tolist()
        return parce_list
    
    def parcellation_100(self):
        d = os.path.join(self.ds_external)
        separator      = ''
        parce_frame    = pd.read_csv(separator.join([d,'/Schaefer2018_100Parcels_7Networks_tab.txt']), delimiter="\t",header=None)
        parce_list     = parce_frame[1].tolist()
        return parce_list
    
    def network_400(self, name_network):
        
        d              = os.path.join(self.ds_external)
        separator      = ''
        parce_frame    = pd.read_csv(separator.join([d,'/Schaefer2018_400Parcels_7Networks_tab.txt']), delimiter=",")
        parce_list     = parce_frame["Var1"].tolist()
        
        mynetlist      = [re.search(name_network,parcel) for parcel in parce_list]
        mynetlocation  = np.array([el is not None for el in mynetlist])
        net_rsn        = np.where(mynetlocation==True)
        net_label      = [parce_list[net_idx] for net_idx in net_rsn[0]]
        
        mynetlist      = [re.search(separator.join(['LH_',name_network]),parcel) for parcel in parce_list]
        mynetlocation  = np.array([el is not None for el in mynetlist])
        left_rsn       = np.where(mynetlocation==True)
        left_label     = [parce_list[net_idx] for net_idx in net_rsn[0]]
        
        mynetlist      = [re.search(separator.join(['RH_',name_network]),parcel) for parcel in parce_list]
        mynetlocation  = np.array([el is not None for el in mynetlist])
        right_rsn      = np.where(mynetlocation==True)
        right_label    = [parce_list[net_idx] for net_idx in net_rsn[0]]
        
        return net_rsn[0],net_label,left_rsn[0],left_label,right_rsn[0],right_label
    
    def network_100(self, name_network):
        
        d              = os.path.join(self.ds_external)
        separator      = ''
        parce_frame    = pd.read_csv(separator.join([d,'/Schaefer2018_100Parcels_7Networks_tab.txt']), delimiter="\t",header=None)
        parce_list     = parce_frame[1].tolist()
        
        mynetlist      = [re.search(name_network,parcel) for parcel in parce_list]
        mynetlocation  = np.array([el is not None for el in mynetlist])
        net_rsn        = np.where(mynetlocation==True)
        net_label      = [parce_list[net_idx] for net_idx in net_rsn[0]]
        
        mynetlist      = [re.search(separator.join(['LH_',name_network]),parcel) for parcel in parce_list]
        mynetlocation  = np.array([el is not None for el in mynetlist])
        left_rsn       = np.where(mynetlocation==True)
        left_label     = [parce_list[net_idx] for net_idx in net_rsn[0]]
        
        mynetlist      = [re.search(separator.join(['RH_',name_network]),parcel) for parcel in parce_list]
        mynetlocation  = np.array([el is not None for el in mynetlist])
        right_rsn      = np.where(mynetlocation==True)
        right_label    = [parce_list[net_idx] for net_idx in net_rsn[0]]
        
        return net_rsn[0],net_label,left_rsn[0],left_label,right_rsn[0],right_label
    
    def cognitive_score(self):
        
        d              = os.path.join(self.ds_external)
        separator      = ''
        score_frame    = pd.read_csv(separator.join([d,'/1000BD_cognitivescore_1visit.csv']),delimiter=";",decimal = ',')
        
        subj_age_CS       = score_frame["Age"].tolist()
        gender_CS         = score_frame["Sex"].tolist()
        subj_ID_CS        = score_frame["ID"].tolist()
        visit_CS          = score_frame["Visit"].tolist()
        processing_speed  = score_frame["TMT_ARW(ProcessingSpeed)"].tolist()
        concept_shifting  = score_frame["TMT_BA(ConceptShifting)"].tolist()
        working_memory    = score_frame["LPS_RRW(ProblemSolving)"].tolist()
        PhonematicFluency = score_frame["RWT_PB2RW(PhonematicFluency)"].tolist()
        SemanticFluency   = score_frame["RWT_SB2RW(SemanticFluency)"].tolist()
        PhonematicFluency_Switch  = score_frame["RWT_PGR2RW(PhonematicFluency_Switch)"].tolist()
        SemanticFluency_Switch    = score_frame["RWT_SSF2RW(SemanticFluency_Switch)"].tolist()
        Vocabulary        = score_frame["AWST03P(Vocabulary)"].tolist()
        
        return subj_age_CS,gender_CS,subj_ID_CS,visit_CS,processing_speed,concept_shifting,working_memory,PhonematicFluency,SemanticFluency,PhonematicFluency_Switch,SemanticFluency_Switch,Vocabulary

    def list_subjects(self):
        d = os.path.join(self.ds_root)
        return [subj for subj in os.listdir(d) if not subj.startswith(".")]
    
    def load_subject_sc(self,subj):
        
        d = os.path.join(self.ds_root)
        separator   = ''
        file_julich = separator.join([d,'/',subj,'/ses-1/SC/',subj,'_SC_Schaefer7NW400p.txt']) 
        file_julich_no_log = separator.join([d,'/',subj,'/ses-1/SC/',subj,'_SC_Schaefer7NW400p_nolog10.txt']) 
        SC = np.loadtxt(file_julich)
        SC_nolog = np.loadtxt(file_julich_no_log)
        return SC, SC_nolog
        
        
    def load_subject_fc(self,subj):
        
        d = os.path.join(self.ds_root)
        separator   = ''
        file_julich = separator.join([d,'/',subj,'/ses-1/FC/',subj,'_meanTS_Schaefer7NW400p.txt'])
        file_julich_noaroma = separator.join([d,'/',subj,'/ses-1/FC/meanTS_RSsw',subj,'_ses-1_bptf_Schaefer7NW400parcels.txt'])
        bold = np.loadtxt(file_julich)
        try:
            bold_noaroma = np.loadtxt(file_julich_noaroma)
            err_message  = 0
        except OSError:
            print('OS-ERROR')
            bold_noaroma = bold
            err_message  = 1
        
        return bold, bold_noaroma, err_message
    
    def load_subject_fc_100(self,subj):
        
        d = os.path.join(self.ds_root)
        
        separator   = ''
        file_julich = separator.join([d,'/',subj,'/ses-1/FC/',subj,'_meanTS_GS_bptf_Schaefer100_7NW.txt'])
        bold        = np.loadtxt(file_julich)
        
        return bold
    
    def load_subject_sc_100(self,subj):
        
        d = os.path.join(self.ds_root)
        
        separator          = ''
        file_julich        = separator.join([d,'/',subj,'/ses-1/SC/',subj,'_SC_Schaefer_7NW100p.txt']) 
        file_julich_no_log = separator.join([d,'/',subj,'/ses-1/SC/',subj,'_SC_Schaefer7NW100p_nolog10.txt']) 
        SC                 = np.loadtxt(file_julich)
        SC_nolog           = np.loadtxt(file_julich_no_log)
        
        return SC, SC_nolog
    
    def load_subject_fc_ses2(self,subj):

        d = os.path.join(self.ds_root)
        separator   = ''
        
        file_fc_400 = separator.join([d,'/',subj,'/ses-2/FC/',subj,'_meanTS_GS_bptf_Schaefer400_7NW.txt'])
        file_fc_100 = separator.join([d,'/',subj,'/ses-2/FC/',subj,'_meanTS_GS_bptf_Schaefer100_7NW.txt']) 
        
        try:
            bold_400     = np.loadtxt(file_fc_400)
            bold_100     = np.loadtxt(file_fc_100)
            err_message  = 0
        except OSError:
            print('NO SESSION 2')
            bold_400     = []
            bold_100     = []
            err_message  = 1

        return bold_400, bold_100, err_message

    def metadata(self):
        d            = os.path.join(self.ds_external)
        d_subj       = os.path.join(self.ds_root)
        
        separator    = ''
        patient      = pd.read_csv(separator.join([d,'/HBP_Descriptives.csv']), delimiter=";",decimal = ',')
        subj_age     = patient["Age"].tolist()
        gender       = patient["Sex"].tolist()
        education    = patient["Education"].tolist()
        subj_ID      = patient["ID_HBP"].tolist()
        
        subjs        = [subj for subj in os.listdir(d_subj) if not subj.startswith(".")]
        
        separator      = ''
        patient        = pd.read_csv(separator.join([d,'/HBP_Descriptives.csv']), delimiter=";",decimal = ',')
        patient        = patient.set_index('ID_HBP')
        patient_sort   = patient.loc[subjs]
        age_sort       = patient_sort["Age"].tolist()
        gender_sort    = patient_sort["Sex"].tolist()
        education_sort = patient_sort["Education"].tolist()
    
        return subj_age,gender,education,subj_ID,age_sort,gender_sort,education_sort,subjs
    
    def metadata_cognitive(self):
        
        d            = os.path.join(self.ds_external)
        d_subj       = os.path.join(self.ds_root)
    
        separator      = ''
        patient       = pd.read_csv(separator.join([d,'/1000Brains_CognitivePerformance/1000Brains_CognitivePerformance_2nd_Visit.csv']),delimiter=";",decimal = ',')
        subj_age_2     = patient["Age"].tolist()
        gender_2       = patient["Sex"].tolist()
        subj_ID_2      = patient["ID"].tolist()
        visit_2        = patient["Visit"].tolist()
    
        return subj_age_2,gender_2,subj_ID_2,visit_2
    
    def load_mask_100(self,unix):
        
        d = os.path.join(self.ds_external)
        
        separator          = ''
        file_mask          = separator.join([d,'/mask/mask',
                                      '_idx',f"{unix}",'_100.npz']) 
        mask_data          = np.load(file_mask)
        mask               = mask_data["mask"]
        
        return mask
    
class Julich_interface:

    def collect_dataset(self,inpath):

        input_data    = np.load(inpath)
        Bold_data     = input_data['Bold_data']
        Bold_time     = input_data['Bold_time']

        return Bold_data,Bold_time,input_data
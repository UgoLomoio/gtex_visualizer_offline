"""
*************************************************************************************************************
*                                                                                                           *
*   GTex data visualizer developed by Ugo Lomoio at Magna Graecia University of Catanzaro                   *
*                                                                                                           *
*                           Dash App with Server Side callbacks                                             *
*                                                                                                           *
*************************************************************************************************************
"""

import pandas as pd 
import scipy.stats as stats
from make_plots import *
import os 

def single_dd_values_handler(gencode_id, gene_name, tissue, filter):

    if tissue == "All":
        #fig_pie = plot_gene_data(gencode_id, gene_name)     
        print("Extracting tpms data")
        if filter == "NoFilters":

            fig_violin, curr_data = plot_by_gene(gencode_id, gene_name)
            return curr_data

        elif filter == "DivideByGender":
                        
            fig_violin, curr_data = plot_by_gene_and_gender(gencode_id, gene_name) 
            return curr_data

        else: 

                 
            return None

    else:
        
        print("Extracting tpms data with filter {}".format(filter))

        if filter == "NoFilters":
                
            fig_violin, curr_data = plot_by_gene_and_tissue(gencode_id, gene_name, tissue)            
            return curr_data

        elif filter == "DivideByGender":

            fig_violin, curr_data = plot_by_gene_and_gender_and_tissue(gencode_id, gene_name, tissue)
            return curr_data

        elif filter == "DivideByAge":

            fig_violin, curr_data = plot_by_gene_and_tissue_and_age(gencode_id, gene_name, tissue)
            return curr_data

        else: #"Divide by Gender and Age"

            fig_violin, curr_data = plot_by_gene_tissue_age_and_gender(gencode_id, gene_name, tissue)
            return curr_data

if __name__ == "__main__":

    print("Start setup")
    all_methods = ["None", "with_labels", "betweenness_centrality", "closeness_centrality", "degree_centrality", "eigenvector_centrality", "community_louvain", "community_leiden", "spectral_clustering"] # "community_girvan_newmann"
    all_tissues = ['All', 'Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland',
     'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Bladder',
     'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
     'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
     'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9',
     'Brain_Hippocampus', 'Brain_Hypothalamus',
     'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia',
     'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra',
     'Breast_Mammary_Tissue', 'Cells_Cultured_fibroblasts',
     'Cells_EBV-transformed_lymphocytes', 'Cervix_Ectocervix',
     'Cervix_Endocervix', 'Colon_Sigmoid', 'Colon_Transverse',
     'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa',
     'Esophagus_Muscularis', 'Fallopian_Tube', 'Heart_Atrial_Appendage',
     'Heart_Left_Ventricle', 'Kidney_Cortex', 'Kidney_Medulla', 'Liver', 'Lung',
     'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary',
     'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic',
     'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen',
     'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood']
    wanted_tissues = ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum',
     'Artery_Aorta', 'Artery_Coronary', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
     'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
     'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9',
     'Brain_Hippocampus', 'Brain_Hypothalamus',
     'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia',
     'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra',
     'Heart_Atrial_Appendage',
     'Heart_Left_Ventricle', 'Kidney_Cortex', 'Kidney_Medulla', 'Liver', 'Lung',
     'Spleen']
    wanted_genes = {}
    with open('ListadiGeni.tex', 'r') as f:
    
        lines = f.readlines()
        for line in lines: 
            wanted_genes[line.split()[1]] = line

    with open("ENSG_to_ENSP.txt", "r") as f:
        ensembl_gene_protein_mapping = dict(eval(f.read()))

    
    with open("all_genes_dict.txt", "r") as f:
        all_genes = dict(eval(f.read()))

    print("Start Plots")
    #all_filters = ["NoFilters", "DivideByGender", "DivideByAge", "DivideByGenderAndAge"]
    all_filters = ["DivideByGenderAndAge"]
    all_ages = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]
    all_genders = ["male", "female"]
    bg = "white" #'rgb(17, 17, 17)'
    txt_color = "black"
    genes_not_in_gtex = []
    for gene_name in wanted_genes.keys():
        print(gene_name)
        if not os.path.exists("TPMs/{}".format(gene_name)):
            os.makedirs("TPMs/{}".format(gene_name))
        files_done = os.listdir("TPMs/{}".format(gene_name))
        if gene_name in all_genes.keys():
            print("exists")
            gencode_id = all_genes[gene_name]
            print(gencode_id)
            for tissue in wanted_tissues:
                print(tissue)
                done_pie = False 
                for filter in all_filters:
                    tmp_filename = "TPMs/{}/tpms_Gene{}_Tissue{}_{}.csv".format(gene_name, gene_name, tissue.replace("_", ""), filter)
                    if tmp_filename.split("/")[2] not in files_done:
                        print(tmp_filename+" not done")
                        df = single_dd_values_handler(gencode_id, gene_name, tissue, filter)
                        if df is not None:
                            #fig_violin.write_image("Plots/Violinplot_Gene{}_Tissue{}_Filter{}.png".format(gene_name, tissue, filter))
                            df.to_csv(tmp_filename)
                    else: 
                        print(tmp_filename+" done")
        else:
            print("doesn't exists")
            genes_not_in_gtex.append(gene_name)

    with(open("genes_not_in_getx_db.txt", "w")) as f:
        for gene_name in genes_not_in_gtex:
            f.write(gene_name)

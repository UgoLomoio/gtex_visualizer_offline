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
from make_plots import *
#import methods
import scipy.stats as stats
import os 

def single_dd_values_handler(gencode_id, gene_name, tissue, filter):

    table_title = "Anova, Kruskal and Shapiro analysis results for: Gene '{}', Tissue '{}' and Filter '{}'".format(gene_name, tissue, filter)

    if tissue == "All":
        fig_pie = plot_gene_data(gencode_id, gene_name)     
        print("Creating violin plots")
        if filter == "NoFilters":

            fig_violin, curr_data = plot_by_gene(gencode_id, gene_name)
                
            ys = [fig_violin.data[i]["y"] for i in range(len(fig_violin.data))]
            ys_shapiro = [elem for y in ys for elem in y]
            shapiro_test = stats.shapiro(ys_shapiro)
            fvalue_shapiro = shapiro_test.statistic
            pvalue_shapiro = shapiro_test.pvalue

        
            fvalue_anova, pvalue_anova = stats.f_oneway(ys[0], ys[1], ys[2], ys[3], ys[4], ys[5], ys[6], ys[7], ys[8], ys[9], ys[10],
                                                ys[11], ys[12], ys[13], ys[14], ys[15], ys[16], ys[17], ys[18], ys[19], ys[20],
                                                ys[21], ys[22], ys[23], ys[24], ys[25], ys[26], ys[27], ys[28], ys[29], ys[30],
                                                ys[31], ys[32], ys[33], ys[34], ys[35], ys[36], ys[37], ys[38], ys[39], ys[40],
                                                ys[41], ys[42], ys[43], ys[44], ys[45], ys[46], ys[47], ys[48], ys[49], ys[50],
                                                ys[51], ys[52], ys[53])
            fvalue_kruskal, pvalue_kruskal = stats.kruskal(ys[0], ys[1], ys[2], ys[3], ys[4], ys[5], ys[6], ys[7], ys[8], ys[9], ys[10],
                                                ys[11], ys[12], ys[13], ys[14], ys[15], ys[16], ys[17], ys[18], ys[19], ys[20],
                                                ys[21], ys[22], ys[23], ys[24], ys[25], ys[26], ys[27], ys[28], ys[29], ys[30],
                                                ys[31], ys[32], ys[33], ys[34], ys[35], ys[36], ys[37], ys[38], ys[39], ys[40],
                                                ys[41], ys[42], ys[43], ys[44], ys[45], ys[46], ys[47], ys[48], ys[49], ys[50],
                                                ys[51], ys[52], ys[53])
            df = pd.DataFrame([["Anova", fvalue_anova, pvalue_anova], ["Kruskal", fvalue_kruskal, pvalue_kruskal], ["Shapiro", fvalue_shapiro, pvalue_shapiro]], columns = ["", "f_value", "p_value"])
            dict_data = df.to_dict('rows')
                        
           
            return fig_violin, fig_pie, table_title, dict_data

        elif filter == "DivideByGender":
                        
            fig_violin, curr_data = plot_by_gene_and_gender(gencode_id, gene_name)
                
            ys = [fig_violin.data[i]["y"] for i in range(len(fig_violin.data))]
            ys_shapiro = [elem for y in ys for elem in y]
            shapiro_test = stats.shapiro(ys_shapiro)
            fvalue_shapiro = shapiro_test.statistic
            pvalue_shapiro = shapiro_test.pvalue

            fvalue_anova, pvalue_anova = stats.f_oneway(ys[0], ys[1], ys[2], ys[3], ys[4], ys[5], ys[6], ys[7], ys[8], ys[9], ys[10],
                                                ys[11], ys[12], ys[13], ys[14], ys[15], ys[16], ys[17], ys[18], ys[19], ys[20],
                                                ys[21], ys[22], ys[23], ys[24], ys[25], ys[26], ys[27], ys[28], ys[29], ys[30],
                                                ys[31], ys[32], ys[33], ys[34], ys[35], ys[36], ys[37], ys[38], ys[39], ys[40],
                                                ys[41], ys[42], ys[43], ys[44], ys[45], ys[46], ys[47], ys[48], ys[49], ys[50],
                                                ys[51], ys[52], ys[53], ys[54], ys[55], ys[56], ys[57], ys[58], ys[59], ys[60],
                                                ys[61], ys[62], ys[63], ys[64], ys[65], ys[66], ys[67], ys[68], ys[69], ys[70],
                                                ys[71], ys[72], ys[73], ys[74], ys[75], ys[76], ys[77], ys[78], ys[79], ys[80],
                                                ys[81], ys[82], ys[83], ys[84], ys[85], ys[86], ys[87], ys[88], ys[89], ys[90],
                                                ys[91], ys[92], ys[93], ys[94], ys[95], ys[96], ys[97], ys[98], ys[99], ys[100],
                                                ys[101], ys[102], ys[103], ys[104], ys[105], ys[106], ys[107])
            fvalue_kruskal, pvalue_kruskal = stats.kruskal(ys[0], ys[1], ys[2], ys[3], ys[4], ys[5], ys[6], ys[7], ys[8], ys[9], ys[10],
                                                ys[11], ys[12], ys[13], ys[14], ys[15], ys[16], ys[17], ys[18], ys[19], ys[20],
                                                ys[21], ys[22], ys[23], ys[24], ys[25], ys[26], ys[27], ys[28], ys[29], ys[30],
                                                ys[31], ys[32], ys[33], ys[34], ys[35], ys[36], ys[37], ys[38], ys[39], ys[40],
                                                ys[41], ys[42], ys[43], ys[44], ys[45], ys[46], ys[47], ys[48], ys[49], ys[50],
                                                ys[51], ys[52], ys[53], ys[54], ys[55], ys[56], ys[57], ys[58], ys[59], ys[60],
                                                ys[61], ys[62], ys[63], ys[64], ys[65], ys[66], ys[67], ys[68], ys[69], ys[70],
                                                ys[71], ys[72], ys[73], ys[74], ys[75], ys[76], ys[77], ys[78], ys[79], ys[80],
                                                ys[81], ys[82], ys[83], ys[84], ys[85], ys[86], ys[87], ys[88], ys[89], ys[90],
                                                ys[91], ys[92], ys[93], ys[94], ys[95], ys[96], ys[97], ys[98], ys[99], ys[100],
                                                ys[101], ys[102], ys[103], ys[104], ys[105], ys[106], ys[107])
                
            df = pd.DataFrame([["Anova", fvalue_anova, pvalue_anova], ["Kruskal", fvalue_kruskal, pvalue_kruskal], ["Shapiro", fvalue_shapiro, pvalue_shapiro]], columns = ["", "f_value", "p_value"])
            dict_data = df.to_dict('rows')
                        
           
            return fig_violin, fig_pie, table_title, dict_data

        else:

                        
            return None, None, None, None

    else:
        
        fig_pie = plot_gene_tissue_data(gencode_id, gene_name, tissue) 
        print("Creating violin plots")

        if filter == "NoFilters":
                
            fig_violin, curr_data = plot_by_gene_and_tissue(gencode_id, gene_name, tissue)
            ys = [fig_violin.data[i]["y"] for i in range(len(fig_violin.data))]
            ys_shapiro = [elem for y in ys for elem in y]
            shapiro_test = stats.shapiro(ys_shapiro)
            fvalue_shapiro = shapiro_test.statistic
            pvalue_shapiro = shapiro_test.pvalue
            df = pd.DataFrame([["Anova", "None", "None"], ["Kruskal", "None", "None"], ["Shapiro", fvalue_shapiro, pvalue_shapiro]], columns = ["", "f_value", "p_value"])
            dict_data = df.to_dict('rows')            
            return fig_violin, fig_pie, table_title, dict_data

        elif filter == "DivideByGender":

            fig_violin, curr_data = plot_by_gene_and_gender_and_tissue(gencode_id, gene_name, tissue)
            ys = [fig_violin.data[i]["y"] for i in range(len(fig_violin.data))]
            fvalue_anova, pvalue_anova = stats.f_oneway(ys[0], ys[1])
            fvalue_kruskal, pvalue_kruskal = stats.kruskal(ys[0], ys[1])
            ys_shapiro = [elem for y in ys for elem in y]
            shapiro_test = stats.shapiro(ys_shapiro)
            fvalue_shapiro = shapiro_test.statistic
            pvalue_shapiro = shapiro_test.pvalue
            df = pd.DataFrame([["Anova", fvalue_anova, pvalue_anova], ["Kruskal", fvalue_kruskal, pvalue_kruskal], ["Shapiro", fvalue_shapiro, pvalue_shapiro]], columns = ["", "f_value", "p_value"])
            dict_data = df.to_dict('rows')
            return fig_violin, fig_pie, table_title, dict_data

        elif filter == "DivideByAge":

            fig_violin, curr_data = plot_by_gene_and_tissue_and_age(gencode_id, gene_name, tissue)
            ys = [fig_violin.data[i]["y"] for i in range(len(fig_violin.data))]
            fvalue_anova, pvalue_anova = stats.f_oneway(ys[0], ys[1], ys[2], ys[3], ys[4])
            fvalue_kruskal, pvalue_kruskal = stats.kruskal(ys[0], ys[1], ys[2], ys[3], ys[4])
            ys_shapiro = [elem for y in ys for elem in y]
            shapiro_test = stats.shapiro(ys_shapiro)
            fvalue_shapiro = shapiro_test.statistic
            pvalue_shapiro = shapiro_test.pvalue
            df = pd.DataFrame([["Anova", fvalue_anova, pvalue_anova], ["Kruskal", fvalue_kruskal, pvalue_kruskal], ["Shapiro", fvalue_shapiro, pvalue_shapiro]], columns = ["", "f_value", "p_value"])
            dict_data = df.to_dict('rows')
            return fig_violin, fig_pie, table_title, dict_data

        else: #"Divide by Gender and Age"

            fig_violin, curr_data = plot_by_gene_tissue_age_and_gender(gencode_id, gene_name, tissue)
            ys = [fig_violin.data[i]["y"] for i in range(len(fig_violin.data))]
            fvalue_anova, pvalue_anova = stats.f_oneway(ys[0], ys[1], ys[2], ys[3], ys[4], ys[5], ys[6], ys[7], ys[8], ys[9])
            fvalue_kruskal, pvalue_kruskal = stats.kruskal(ys[0], ys[1], ys[2], ys[3], ys[4], ys[5], ys[6], ys[7], ys[8], ys[9])
            ys_shapiro = [elem for y in ys for elem in y]
            shapiro_test = stats.shapiro(ys_shapiro)
            fvalue_shapiro = shapiro_test.statistic
            pvalue_shapiro = shapiro_test.pvalue
            df = pd.DataFrame([["Anova", fvalue_anova, pvalue_anova], ["Kruskal", fvalue_kruskal, pvalue_kruskal], ["Shapiro", fvalue_shapiro, pvalue_shapiro]], columns = ["", "f_value", "p_value"])
            dict_data = df.to_dict('rows')
            return fig_violin, fig_pie, table_title, dict_data

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
    all_filters = ["NoFilters", "DivideByGenderAndAge"]
    all_ages = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]
    all_genders = ["male", "female"]
    bg = "white" #'rgb(17, 17, 17)'
    txt_color = "black"
    genes_not_in_gtex = []
    for gene_name in wanted_genes.keys():
        print(gene_name)
        if not os.path.exists("Plots/{}".format(gene_name)):
            os.makedirs("Plots/{}".format(gene_name))
        files_done = os.listdir(os.getcwd()+"/Plots/"+gene_name)
        if gene_name in all_genes.keys():
            print("exists")
            gencode_id = all_genes[gene_name]
            print(gencode_id)
            for tissue in all_tissues:
                print(tissue)
                done_pie = False 
                for filter in all_filters:
                    violin_filename = "Plots/{}/Violinplot_Gene{}_Tissue{}_Filter{}.html".format(gene_name, gene_name, tissue, filter)
                    pie_filename = "Plots/{}/Pieplot_Gene{}_Tissue{}.html".format(gene_name, gene_name, tissue, filter)
                    if violin_filename.split("/")[2] not in files_done:
                        print(violin_filename+" not done")
                        fig_violin, fig_pie, table_title, dict_data = single_dd_values_handler(gencode_id, gene_name, tissue, filter)
                        if fig_violin is not None:
                            #fig_violin.write_image("Plots/Violinplot_Gene{}_Tissue{}_Filter{}.png".format(gene_name, tissue, filter))
                            fig_violin.write_html(violin_filename)
                        if fig_pie is not None:
                            if pie_filename.split("/")[2] not in files_done:
                                if not done_pie:
                                    print(pie_filename+" not done")
                                    #fig_pie.write_image("Plots/Pieplot_Gene{}_Tissue{}.png".format(gene_name, tissue, filter))
                                    fig_pie.write_html(pie_filename)
                                    done_pie = True
                                    print("done")
                            else:
                                print(pie_filename+" done")
                    else: 
                        print(violin_filename+" done")
        else:
            print("doesn't exists")
            genes_not_in_gtex.append(gene_name)

    with(open("genes_not_in_getx_db.txt", "w")) as f:
        for gene_name in genes_not_in_gtex:
            f.write(gene_name)



# 307genes × 60tissues × 5mb = 92000MB = 100GB of plots
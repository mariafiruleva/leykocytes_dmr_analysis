rule data_vizualisation:
    input: "all_cpgs.10x_leukocytes_everywhere_n=36.txt"
    output: "violin_plot.png"
    shell: "python3.6 data_vizualisation.py"


rule data_preparation:
    input: "all_cpgs.10x_leukocytes_n=41.csv"
    output: "all_cpgs.10x_leukocytes_everywhere_n=36.txt"
    shell: "python3.6 data_preparation.py"

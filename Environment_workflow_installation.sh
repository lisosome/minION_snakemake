#!/usr/bin/env bash
set -e
# Commands to install miniconda and/or snakemake 6.15.5 for pipelines or both


install_miniconda=false
install_snakemake=false
install_all=false

printHelp(){
    echo "Environment_workflow_installation.sh"
    echo 
    echo "Usage: Environment_workflow_installation.sh [options]"
    echo "Options:"
    echo "-h, --help - Show this help."
    echo "-m, --miniconda - Install only miniconda3"
    echo "-s, --snakemake - Install only snakemake"
    echo "-a, --all - Install both miniconda and snakemake"
}

if [[ $# -eq 0 ]];then 
    printHelp
    exit 0
fi

if [[ $# -gt 1 ]];then 
    echo "Only one flag accepted. Please select a proper flag"
    printHelp
    exit 1
fi



inst_conda(){
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh &&
    ~/miniconda3/bin/conda init bash && source ~/.bashrc
}

inst_snake(){
    conda create -y -n snakemake -c conda-forge python=3.10.2 &&
    source activate snakemake && conda install -y -c bioconda -c conda-forge snakemake &&
    conda install -y -c bioconda -c conda-forge --force-reinstall snakemake=6.15.5 &&
    conda install -y -c conda-forge --force-reinstall tabulate=0.8
}

while [[ $# -gt 0 ]];do
    arg=$1
    case ${arg} in
        
        -h|--help)
        printHelp
        exit 0
        ;;
        
        -m|--miniconda)
        install_miniconda=true
        break
        ;;

        -s|--snakemake)
        install_snakemake=true
        break
        ;;

        -a|--all)
        install_all=true
        break
        ;;
        *)
        echo "Unknown option: $1"
        printHelp
        exit 1
    esac
    shift
done


#echo "Install Miniconda: $install_miniconda"
#echo "Install Snakemake: $install_snakemake"
#echo "Install all: $install_all"

if [[ $install_all == true && ($install_miniconda == true || $install_snakemake == true) ]];then
    script=$(realpath $0)
    echo "-a | --all flag must be the only argument. Please run the following command: ${script} -a"
    exit 1
elif [[ $install_miniconda == true && $install_all == false && $install_snakemake == false ]];then
        #echo "this choice install_miniconda=true && install_all=false && install_snakemake=false"
        inst_conda
elif [[  $install_miniconda == false && $install_all == true && $install_snakemake == false ]];then
        inst_conda && inst_snake
        #echo "this choice install_miniconda=false && install_all=true && install_snakemake=false"
elif 
    [[ $install_miniconda == false && $install_all == false && $install_snakemake == true ]];then
        inst_snake
        #echo "this choice install_miniconda=false && install_all=false && install_snakemake=true"
else
    echo "Select a valid flag"
    printHelp
    exit 0
fi












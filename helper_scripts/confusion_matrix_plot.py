import sys
sys.path.insert(0,r'..') # to allow import packages in parent folder
import pandas as pd
import run_grep
import glob
import matplotlib.pyplot as plt
import numpy as np
import random

# helper method to create downsampled raw reads and index files
def prep_downsampling_ratio_data(raw_reads_folder, ratios):
    fastq_R1_file = rf"{raw_reads_folder}/Undetermined_S0_L001_R1_001.trim.fastq"
    fastq_R2_file = rf"{raw_reads_folder}/Undetermined_S0_L001_R2_001.trim.fastq"

    fastq_I1_file = rf"{raw_reads_folder}/Undetermined_S0_L001_I1_001.fastq"
    fastq_I2_file = rf"{raw_reads_folder}/Undetermined_S0_L001_I2_001.fastq"

    for ratio in ratios:
        command = (rf"reformat.sh in1={fastq_R1_file} in2={fastq_R2_file} out1={raw_reads_folder}/out{ratio}.R1.fastq out2={raw_reads_folder}/out{ratio}.R2.fastq "
                    f"samplerate={ratio} sampleseed=77")
        run_grep.exec_grep_shell(command)
        command = (rf"reformat.sh in1={fastq_I1_file} in2={fastq_I2_file} out1={raw_reads_folder}/out{ratio}.I1.fastq out2={raw_reads_folder}/out{ratio}.I2.fastq "
                    f"samplerate={ratio} sampleseed=77")
        run_grep.exec_grep_shell(command)

# helper method to set up downsampling folders (creat the folder and copy 2 files into it)
def set_up_directory(qc_output_folder, ratios):
    for ratio in ratios:
        command = (f"mkdir -p {qc_output_folder}/output{ratio} && "
            f"cp {qc_output_folder}/settings.ini {qc_output_folder}/Undetermined.paired.files {qc_output_folder}/output{ratio}/")
        run_grep.exec_grep_shell(command)

# helper method to create batch files for each downsampling folder
def write_batch_file(qc_output_folder, ratios):
    for ratio in ratios:
        batch_file = f"out{ratio}.R1.fastq\tout{ratio}.R2.fastq\tout{ratio}.I1.fastq\tout{ratio}.I2.fastq"
        command = (f"cd {qc_output_folder}/output{ratio} && "
            f"echo {batch_file} > Undetermined.paired.files")
        run_grep.exec_grep_shell(command)

# convinient method to set up and call all required methods
def prep_downsampling():
	qc_output_folder = r"/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_QC_pipeline/M3235_22_019"
	ratios = [0.5,0.2,0.1,0.05,0.02,0.01,0.005]
	raw_reads_folder = r"/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/cgMLST_pilot/M3235-22-019_fastq"
	#1
	prep_downsampling_ratio_data(raw_reads_folder,ratios)
	#2
	set_up_directory(qc_output_folder,ratios)
	#3
	write_batch_file(qc_output_folder,ratios)
	#4
	## you will still need to edit the 'output' section of settings.ini in each folder
	#5 
	## kick off the qc pipeline for each downsampling data set


def make_df_for_plot_confusion_matrix(dir, metric):
    '''
    this method check all the files in the given directory. Read each file into a DataFrame, and make a new
    DataFrame by concating the intended 'metric' column of each DataFrame

    Parameters
    ----------
    dir: String (folder name for the input files)
    metric: String (the intended column name)

    Returns
    ----------
    a DataFrame (concated with all the 'metric' columns in the input files)

    '''
    df = None
    for file in sorted(glob.glob(rf'{dir}/*.txt')):
        if df is None:
            df = pd.read_csv(file, sep='\t', index_col=[0])[metric]
        else:
            new_df = pd.read_csv(file, sep='\t', index_col=[0])[metric]
            df = pd.concat([df, new_df], axis=1)

    return df



def plot_confusion_matrix(df, metric, sample_list, x_axis, x, title, saved_file):
    '''
    this method plots one specific given metric (from confusion matrix) for the downsampling data set.

    Parameters
    ----------
    df: DataFrame (row: samples, column: the specific metric value from different downsampling data)
    metric: String (metric from confusion_matrix. Ex. sensitivity (TP/P))
    sample_list: list (list of the samples we're studying)
    x_axis: list (x axis value)
    x: list (intended x axis value. for xticks())
    title: String (tile of the plot)
    saved_file: String (name of the file to be saved)

    '''
    # x_axis = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    # # x_axis = np.log10([x+0.1 for x in x_axis])
    # x = ['(0.1)','(0.2)','(0.3)','(0.4)','(0.5)','(0.6)','(0.7)','(0.8)','(0.9)','original(1)']
    color_list = ['b','g','k','r','m','y']
    marker_list = ['o','+',',','^','.','<']

    control_list = ['Blank_IRP1','Blank_IRP2','Water1']
    sample_list = [i for i in sample_list if i not in control_list]

    fig,ax = plt.subplots()
    for idx, sample in enumerate(random.sample(sample_list,5)):
        ax.plot(x_axis, df.loc[sample].values, color=color_list[idx], marker=marker_list[idx], label=sample)
    ax.plot(x_axis, df.mean().values, color=color_list[idx+1], marker=marker_list[idx+1], label='mean')

    ax.set_xlabel("subsampling ratio (log10)",fontsize=10)
    ax.set_ylabel(metric, color="red",fontsize=9)
    plt.xticks(x_axis, x, fontsize=7)
    plt.legend(loc=7)
    plt.title(title)
    plt.show()
    fig.savefig(saved_file, dpi=1600)

# convenient method to set up and make the cm plot
def plot_cm():
	x_axis = [0.02,0.04,0.06,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
	x_axis = np.log10(x_axis)
	x = ['(0.02)','(0.04)','(0.06)','(0.08)',
		'(0.1)','(0.2)','(0.3)','(0.4)','(0.5)','(0.6)','(0.7)','(0.8)','(0.9)','(1)']
	metric = 'sensitivity (TP/P)'
	# folder of all the confusion_matrix txt files
	cm_dir = r'/scicomp/home-pure/qtl7/HMAS-QC-Pipeline/helper_scripts/confusion_matrix_downsampling2094660'
	title = f'Based on initial reads count 523,650'
	saved_file = f'confusion_matrix_downsampling2094660.pdf'
	plot_confusion_matrix(make_df_for_plot_confusion_matrix(cm_dir, metric), metric, sample_list, x_axis, x, title, saved_file)


def plot_confusion_matrix_ROC(df_x, df_y, title, saved_file):
    '''
    this method plots ROC plot for the downsampling data set.

    Parameters
    ----------
    df_x: DataFrame (row: samples, column: true positive rate value from different downsampling data)
    df_y: DataFrame (row: samples, column: false positive rate value from different downsampling data)
    title: String (tile of the plot)
    saved_file: String (name of the file to be saved)

    '''
    color_list=['b','g','y']
    marker_list=['^','+','<']
    # pick these 3 samples as representative samples
    picked_list = ['08_0810','2010K_1649','2013K_1828']

    fig,ax = plt.subplots()
    for idx, sample in enumerate(picked_list):
        x_axis = list(df_x.loc[sample].values)
        y_axis = list(df_y.loc[sample].values)
        x_axis.insert(0,0)
        y_axis.insert(0,0)
        x_axis.append(1)
        y_axis.append(1)

        ax.plot(x_axis, y_axis, color=color_list[idx], marker=marker_list[idx], label=sample)

    plot_confusion_matrix_mean_ROC(df_x, df_y, title, saved_file, ax, fig)


# this method plots ROC plot for the mean value of downsampling data set.
# it is a special version of plot_confusion_matrix_ROC()
def plot_confusion_matrix_mean_ROC(df_x, df_y, title, saved_file, ax=None, fig=None):

    x_axis = list(df_x.mean(axis=0).values)
    y_axis = list(df_y.mean(axis=0).values)

    # add 0, 1 to the axis
    x_axis.insert(0,0)
    y_axis.insert(0,0)
    x_axis.append(1)
    y_axis.append(1)
    # calculate auc value
    auc = np.trapz(y_axis,x_axis)

    if ax is None:
        fig,ax = plt.subplots()

    ax.plot(x_axis, y_axis,
        color='b', linestyle='dashed', marker='.', markerfacecolor='r', markersize=10,
        label=f'mean \nAUC={auc:.3f}')

    # add annotation to the points
    ### need to customize the list to your own data !
    ratio = [0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,'pool_with_026 (019 data pooled with 026)']
    # ratio = [0.1,0.2,0.3,0.4,0.5,0.6]
    for idx, text in enumerate(ratio, start=1):
        ax.annotate(text, (x_axis[idx], y_axis[idx]),fontsize=8,
             xytext=(x_axis[idx]+0.05*idx, y_axis[idx]-0.3+0.02*idx), 
             arrowprops = dict(arrowstyle="wedge,tail_width=0.2", alpha=0.1))

    # add annotation 'down-sampling ratio'
    idx -= 1 # annotate at the 2nd to last point
    ax.annotate('down-sampling ratio', (x_axis[idx], y_axis[idx]),fontsize=8,
        xytext=(x_axis[idx]+0.05*idx, y_axis[idx]-0.33+0.02*idx))

    ax.set_ylabel("True Positive Rate",fontsize=10)
    ax.set_xlabel("False Positive Rate",fontsize=9)
    plt.legend(loc=7)
    plt.title(title)
    plt.show()
    fig.savefig(saved_file, dpi=1600)

    distance_dict = {(x,y):x**2+(y-1)**2 for x,y in zip(x_axis,y_axis)}
    print (dict(sorted(distance_dict.items(), key=lambda x:x[1])))

    # find the point with min distance to (0,1)
    distance_list = [x**2+(y-1)**2 for x,y in zip(x_axis,y_axis)]
    print (distance_list.index(min(distance_list)))
    print (min(distance_list))
    print (np.argsort(distance_list))


# convenient method to set up and make the cm ROC plot
def plot_cm_ROC():
	# dir = r'/scicomp/home-pure/qtl7/HMAS-QC-Pipeline/helper_scripts/confusion_matrix_downsampling2094660'
	# # dir = r'/scicomp/home-pure/qtl7/HMAS-QC-Pipeline/helper_scripts/confusion_matrix_downsampling364484'
	# title = f'ROC on data with initial reads count 523,650'
	# saved_file = f'ROC_downsampling2094660.pdf'
	# df_y = make_df_for_plot_confusion_matrix(dir, 'sensitivity (TP/P)')
	# df_x = 1 - make_df_for_plot_confusion_matrix(dir, 'specificity (TN/N)')
	# plot_confusion_matrix_ROC(df_x, df_y, title, saved_file)

	# dir = r'/scicomp/home-pure/qtl7/HMAS-QC-Pipeline/helper_scripts/confusion_matrix_downsampling364484'
	# title = f'ROC on data with initial reads count 91,121'
	# saved_file = f'ROC_downsampling364484.pdf'
	# df_y = make_df_for_plot_confusion_matrix(dir, 'sensitivity (TP/P)')
	# df_x = 1 - make_df_for_plot_confusion_matrix(dir, 'specificity (TN/N)')
	# plot_confusion_matrix_ROC(df_x, df_y, title, saved_file)

	dir = r'/scicomp/home-pure/qtl7/HMAS-QC-Pipeline/helper_scripts/confusion_matrix_019'
	title = f'ROC on 019 data with initial mean reads count 826,246'
	saved_file = f'ROC_downsampling019.pdf'
	df_y = make_df_for_plot_confusion_matrix(dir, 'sensitivity (TP/P)')
	df_x = 1 - make_df_for_plot_confusion_matrix(dir, 'specificity (TN/N)')
	# plot_confusion_matrix_mean_ROC(df_x, df_y, title, saved_file)
	plot_confusion_matrix_ROC(df_x, df_y, title, saved_file)


if __name__ == "__main__":

	# prep_downsampling()
	# plot_cm()

	plot_cm_ROC()
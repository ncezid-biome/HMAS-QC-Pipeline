import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from operator import truediv
from primer_length import get_primer_length
# pd.set_option('display.max_columns', None)

file1 = 'final_df_1207.txt'
file2 = 'final_df_1130.txt'
file3 = 'final_df_1206.txt'
file4 = 'final_df_1220.txt'

control_list = ['2013K_0676',
				'2013K_1246',
				'2014K_0979',
				'2016K_0878',
				'Blank_IRP1',
				'Blank_IRP2',
				'Water1']


def file_process_venn(file_name):
	'''
	This method reads the 'final_df' file and returns a processed dataframe and dictionary

	Parameters
    ----------
    file_name: one of the final_df file

    Returns tuple of dataframe and dictionary
    -------
	'''

	df = pd.read_csv(file_name, sep='\t', index_col=0)
	df.drop(df.columns[[0,1,2]], axis=1, inplace=True) #remove 'fail'/'pass<=5'/'pass>5' columns
	df = df.replace(r'^\s*$', np.nan, regex=True) #replace empty cells with nan
	df.fillna(0,inplace=True)
	#dictionary of index:list of primers which has 0 value (failed)
	dic = {idx:list(df.columns[(df==0).loc[idx,]]) for idx in df.index if idx not in control_list}

	return (dic,df)


def file_process_bin(file_name):
	'''
	This method reads the 'final_df' file and returns a processed dataframe with binary value

	Parameters
    ----------
    file_name: one of the final_df file

    Returns a dataframe
    -------
	'''

	df = pd.read_csv(file_name, sep='\t', index_col=0)
	df.drop(df.columns[[0,1,2]], axis=1, inplace=True)
	df = df.replace(r'^\s*$', np.nan, regex=True) #replace empty cells with nan
	df.fillna(0,inplace=True)
	df[df > 0] = 1

	# remove rows whose index is in the control list
	return df.filter(items=[idx for idx in df.index if idx not in control_list], axis=0)


# df_1207 = file_process_bin(file1)
# df_1130 = file_process_bin(file2)
# df_1206 = file_process_bin(file3)
# df_1220 = file_process_bin(file4)

# scatter plot of primer length (x-axis) over # of samples where primer pair succeed (y-axis)
def plot_length_dist(df):
	
	y_axis = list(df.sum(axis=0).values)
	dic = get_primer_length()
	x_axis = [dic[primer] for primer in df.columns ]
	# print (x_axis)
	# print (max(x_axis), min(x_axis)) #440, 180
	print (df.columns[np.argmax(x_axis)])
	# print (y_axis)
	# plt.plot(x_axis,y_axis, 'r.')
	plt.scatter(x_axis,y_axis,s=80, alpha=0.12)
	plt.axis([175, 255, 0, 16]) #confine axis to x[175,255] and y[0,16]
	plt.xlabel('amplicon length (bp)')
	plt.ylabel('number of samples where a primer pair succeeded')
	plt.title('cutadapt trimming M347-21-026')
	plt.savefig('read_length_dist_cutadapt_trimming_026_new.pdf', dpi=1600)
	plt.show()


def plot_length_abundance(df, y_axis, axis_range, out_file, y_label, x_label='amplicon length (bp)', 
							title='cutadapt trimming M347-21-026'):
	
	# df = df.filter(items=[idx for idx in df.index if idx not in control_list], axis=0)
	# y_axis = list(df.mean(axis=0).values)
	# y_axis = list(map(truediv, list(df.std(axis=0).values), list(df.mean(axis=0).values)))
	dic = get_primer_length()
	x_axis = [dic[primer] for primer in df.columns ]
	# print (x_axis)
	# print (max(x_axis), min(x_axis)) #440, 180
	print (df.columns[np.argmax(x_axis)])
	# plt.plot(x_axis,y_axis, 'r.')
	plt.scatter(x_axis,y_axis,s=80, alpha=0.12)
	# plt.axis([175, 255, 0, 400]) #confine axis to x[175,255] and y[0,400]
	plt.axis(axis_range) #confine axis to x[175,255] and y[0,4]
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)
	# plt.savefig('CV_length_dist_nocutadapt_trimming_026.pdf', dpi=1600)
	plt.savefig(out_file, dpi=1600)
	plt.show()



# plot histogram (number of samples where primer pair succeed)
def plot_hist(df, x_label='number of samples where a primer pair succeeded', y_label='count', title='cutadapt trimming M347-21-026'):

	df = df.filter(items=[idx for idx in df.index if idx not in control_list], axis=0)
	arr = plt.hist(list(df.sum(axis=0).values), histtype='stepfilled', bins=[i for i in range(17)])
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)
	for i in range(16): #we have 0..15 values
	    plt.text(arr[1][i],arr[0][i]+20,str(int(arr[0][i])))

	plt.savefig('primer_bin_cutadapt_trimming_026_new.pdf', dpi=1600)
	plt.show()


# plot_hist(df_1206)

def plot_hist_abundance(bin_list, bins, out_file, x_label='number of samples where a primer pair succeeded', 
						y_label='count', title='cutadapt trimming M347-21-026'):

	# arr = plt.hist(list(df.mean(axis=0).values), histtype='stepfilled', bins=[i for i in range(17)])
	# print (df.std())
	# bin_list = list(df.mean(axis=0).values)
	# bin_list = list(map(truediv, list(df.std(axis=0).values), list(df.mean(axis=0).values)))
	# print (min(bin_list), max(bin_list))
	arr = plt.hist(bin_list, histtype='stepfilled', bins=bins)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)
	for i in range(bins): #we have 0..19 values
		if int(arr[0][i]) > 0:
			plt.text(arr[1][i]+10,arr[0][i]+10,str(int(arr[0][i])))
			# plt.text(arr[1][i],arr[0][i]+10,str(int(arr[0][i])), fontsize=6, fontfamily='monospace')

	# plt.savefig('primer_bin_cutadapt_trimming_026.pdf', dpi=1600)
	plt.savefig(out_file, dpi=1600)
	plt.show()


# for generating venn diagrams for failed primers for 2 methods
def make_venn(file1, file2):

	dic_1207,df_1207 = file_process_venn(file1)
	dic_1130,df_1130 = file_process_venn(file2)

	for key in dic_1130:
		print (key, f"{file1} #fail {len(dic_1207[key])}, {file2} #fail {len(dic_1130[key])}")
		only_1207_list = list(set(dic_1207[key]) - set(dic_1130[key]))
		print (f"#of same failed primers: {len(dic_1207[key]) - len(only_1207_list)}")
		print (f"in {file1} only: {len(only_1207_list)}")
		print (f"in {file2} only: {len(list(set(dic_1130[key]) - set(dic_1207[key])))}\n")

		only_1130_list = list(set(dic_1130[key]) - set(dic_1207[key]))
		print(f"{file1} failed primers which appears in {file2} list {len([i for i in only_1130_list if i in list(df_1207.columns)])} times")

		print(f"{file2} failed primers which appears in {file1} list {len([i for i in only_1207_list if i in list(df_1130.columns)])} times\n")

# make_venn(file1,file2)

# plot_hist(df_1220)
# plot_length_dist(df_1220)

def print_matched_per_by(df, cond):
	df.drop(df.columns[[0,1,2]], axis=1, inplace=True)
	df = df.replace(r'^\s*$', np.nan, regex=True) #replace empty cells with nan
	df.fillna(0,inplace=True)
	df = df.filter(items=[idx for idx in df.index if idx not in control_list], axis=0)
	print (df[cond].count().sum() / df[df >= 0].count().sum())



if __name__ == "__main__":
	# df_220104 = file_process_bin(f'final_df_220104_outbase026.txt')
	# plot_hist(df_220104)
	# plot_length_dist(df_220104)

	df = pd.read_csv(f'final_df_220104_outbase026.txt', sep='\t', index_col=0)
	print_matched_per_by(df, df == 1)
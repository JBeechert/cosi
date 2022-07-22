import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_data(filepath,sep=',',header=0,index_col=None,parse_dates=[10],infer_datetime_format=True):
    df = pd.read_csv(filepath)
    df.rename(columns={'Dectector':'Detector'}, inplace=True)
    df.rename(columns={'Background Left':'BackgroundLeft'},inplace=True)
    df.rename(columns={'Background Right':'BackgroundRight'},inplace=True)
    df.rename(columns={'Amplitude Error':'AmplitudeError'},inplace=True)
    df.rename(columns={'Center Error':'CenterError'},inplace=True)
    df.rename(columns={'Width Error':'FWHMError'}, inplace=True)
    df.rename(columns={'Background Left Error':'BackgroundLeftError'}, inplace=True)
    df.rename(columns={'Background Right Error':'BackgroundRightError'}, inplace=True)
    return df

def plot_data(source,parameter,side,dets,filepaths):
    '''
    Plot one parameter exported by the GSE's Auto Fit Panel against strip number.
    This function generates 12 plots, one per detector for either AC or DC.

    Arguments
      - source: Radioactive source name
      - parameter: One of the following: Peak,Amplitude,Center,FWHM,BackgroundLeft,BackgroundRight,AmplitudeError,CenterError,FWHMError,BackgroundLeftError,BackgroundRightError
      - side: AC or DC
      - dets: Detectors to plot, given as [startdetector,enddetector]. To plot only detector 0, give [0,0]. To plot all, give [0,11].
      - filepaths: Up to six file paths to the .csv file of interest, named 'fit_parameters(Source_startdate of data collection).csv'. ex) /home/jacqueline/check_widths/Jan2019_Cs137/'fit_parameters(Cs-137_190114).csv'. Note that this is not the default file name of the file exported by the Auto Fit Panel because that file name corresponds to the date and time the file was generated and not to the data contained therein.

    An example of how to run the code: python plot_params.py Cs-137 FWHM AC [0,11] /home/jacqueline/check_widths/Jan2019_Cs138/'fit_parameters(Cs-137_190114).csv' /home/jacqueline/check_widths/Nov2019_Cs137/'fit_parameters(Cs-137_191113).csv'  
    '''

    fig = plt.figure(figsize=(8,8))
    dfs = {}
    colors = ['red','black','blue','green','magenta','cyan']

    for i in range(len(filepaths)):
        filepath = filepaths[i]
        dfs[filepath] = get_data(filepath)

    color_index = 0

    if (parameter == "WidthError") or (parameter == "Width Error"):
        parameter = "FWHMError"

    if not ("Error" in parameter):
        error = parameter+"Error" 

    for filepath,df in dfs.items():
        data_label = filepath[filepath.find("(")+1:filepath.find(")")]	
        startdet = int(dets[dets.find("[")+1])
        enddet = int(dets[dets.find(",")+1:filepath.find("]")])
         
        for d in range(startdet+1,enddet+2):
            ax = fig.add_subplot(4,3,d)
            mask = ((df['Side'].values == side) & (df['Detector'].values == d-1))
            df_masked =  df[mask]
            
            if ("Error" in parameter) or (parameter == "Peak"):
                ax.plot(df_masked['Strip'].values,df_masked['{}'.format(parameter)].values,'.',c=colors[color_index],label=data_label)
            else:
                ax.errorbar(df_masked['Strip'].values,df_masked['{}'.format(parameter)].values,yerr=df_masked[error].values,marker='.',ls='None',c=colors[color_index],label=data_label)
                
                # Wanted to plot without error bars
                #ax.plot(df_masked['Strip'].values,df_masked['{}'.format(parameter)].values,'.',c=colors[color_index],label=data_label)
                
            ax.axhline(np.mean(df_masked['{}'.format(parameter)]),c=colors[color_index],label='Mean {}'.format(parameter))
            ax.set_title('{}, Detector {}, Side {}'.format(source,d-1,side),fontsize='small')
            ax.set_ylabel('{}'.format(parameter),fontsize='small')
            ax.set_xlabel('Strip #',fontsize='small')
        
        color_index += 1

    plt.legend(bbox_to_anchor=(1.,1.),ncol=1)
    plt.subplots_adjust(left=0.03,bottom=0.05,right=0.69,top=0.97,wspace=0.2,hspace=0.33)
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) > 11:
        print('Exceeded max number of files. Please provide six or fewer files. Do you really want to plot more than six data sets? Too messy!')
        sys.exit(0)

    if sys.argv[1] == "--h" or sys.argv[1] == "--help" or sys.argv[1] == "-h" or sys.argv[1] == "-help":
        print('''
    Plot one parameter exported by the GSE's Auto Fit Panel against strip number.
    This function generates 12 plots, one per detector for either AC or DC.

    Arguments
      - source: Radioactive source name
      - parameter: One of the following: Peak,Amplitude,Center,FWHM,BackgroundLeft,BackgroundRight,AmplitudeError,CenterError,FWHMError,BackgroundLeftError,BackgroundRightError
      - side: AC or DC
      - dets: Detectors to plot, given as [startdetector,enddetector]. To plot only detector 0, give [0,0]. To plot all, give [0,11].
      - filepaths: Up to six file paths to the .csv file of interest, named 'fit_parameters(Source_startdate of data collection).csv'. ex) /home/jacqueline/check_widths/Jan2019_Cs137/'fit_parameters(Cs-137_190114).csv'. Note that this is not the default file name of the file exported by the Auto Fit Panel because that file name corresponds to the date and time the file was generated and not to the data contained therein. 
    ''')
    else:
        source = sys.argv[1]
        parameter = sys.argv[2]
        side = sys.argv[3]
        dets = sys.argv[4]
        filepaths = sys.argv[5:]
        plot_data(source,parameter,side,dets,filepaths)

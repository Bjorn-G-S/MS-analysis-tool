'''
Written by Bjørn Gading Solemsli
date: 14.06.2021
updated: 05.05.2023

Affiliation: iCSI, SMN Catalysis, UiO

'''


from typing import Text
from ipywidgets.widgets.interaction import interactive
from ipywidgets.widgets.widget_float import FloatText
from ipywidgets.widgets.widget_selection import Dropdown
from numpy.core.fromnumeric import size
from numpy.core.function_base import linspace
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.signal as spsign
import ipympl
import matplotlib.colors as colors
import numpy as np
import scipy.integrate as sp# trapz, Simps, cumtrapz, romb
import ipywidgets as widgets

from IPython.display import display
from scipy import sparse
from scipy.sparse.linalg import spsolve


class MS_analysis_MTM():
    
    
    
    
############################################################################################
#                                                                                          #
#                                   Initiating the class.                                  #
#         by doing so, a directory for the TPR files to be analysed will be needed         #
#                                                                                          #
############################################################################################
    
    def __init__(self, *args, **kwargs):
        #calibration in the MS
        self.CO_calib_factor = 0.0000000808    #ppm
        self.CO2_44_calib_factor = 0.00000006329    #ppm
        self.DME_46_calib_factor = 0.00000002444    #ppm
        self.MeOH_29_calib_factor = 0.00000003094    #ppm
        self.MeOH_31_calib_factor = 0.00000001846    #ppm





        Main_Direct_MS = r'C:\Users\bjorngso\OneDrive - Universitetet i Oslo\01 Results\Testing\M2M\Raw data'
        # Main_Direct_MS =r'C:\Users\bjorngso\OneDrive - Universitetet i Oslo\01 Results\Testing\ESRF2023\MS Raw'
        
        together_MS = list()
        for root, dirs, files in os.walk(Main_Direct_MS, topdown=False):
            for name in files:
                if name[-4:] == '.asc':
                    joint = ('{}'.format(name), os.path.join(root, name))
                    together_MS.append(joint)
        
        

        direct_MS= widgets.Dropdown(options=together_MS, description='MS experiment:', disabled=False,layout=widgets.Layout(width='90%'),style = {'description_width': 'initial'})
        
        
        w = widgets.interactive(self.getting_TPD_files,file_MS=direct_MS,)
        display(w)

        self.help()


    def help(self,**kwargs):

        print('''



                   THIS CLASS IS USED FOR ANALYZING EXPERIMENTS (MS AND TEMPERATURE DATA) TAKEN AT THE M2M-rig 
                                              use the following commands in the MS_analysis:



        obj.help()                    -        Will give this list over commands, and what their function is.



        
                                               Methane to Methanol experiments:  



        obj.plot_MS()                 -        Old script for analyzing MS data from a MTM catalytic test
        
        obj.plot_MS_Area()            -        Defining the Area for the MS masses

        obj.Area_to_csv()             -        Used to convert the area peaks to .csv so that it can by analyzed outside the script. 
        
        obj.testing_resutls()         -        When Area is defined form obj.plot_MS_Area, this will do the calcualtion to get the yield, selectivity etc.




        
        ''')
        





    
    
    def getting_TPD_files(self,**kwargs):
        default_MS = r'C:\Users\bjorngso\OneDrive - Universitetet i Oslo\01 Results\Testing\M2M\Raw data\210525_A2-CuZSM5_CH4-TPR.asc'
        
        Filename_MS = kwargs.get('file_MS', default_MS)
        self.Filename_MS = Filename_MS
        
        ### Loading MS files form ASCII type (MTM)
        df_mz = pd.read_csv(Filename_MS, skiprows=7, sep='\t', decimal=",").iloc[:,:-1]
        header = pd.read_csv(Filename_MS, skiprows=5, nrows=0, sep='\t\t\t', engine='python').columns.to_list()
        columns = pd.MultiIndex.from_product([header, ['Time','Time Relative','Ion Current']], names=['m/z', 'data'])
        df_mz.columns = columns


        # df_mz['Ion Current','data'].astype(float)
    
        # Pressure Data
        # df_pressure = pd.read_csv(Filename_MS, skiprows=7, sep='\t', decimal=',').iloc[:, -2]
        # pressure_header = pd.read_csv(Filename_MS, skiprows=5, nrows=0, sep='\t\t\t', engine='python').columns.to_list()[-1:]
        # columns_pressure = pd.MultiIndex.from_product([pressure_header, ['Pressure']], names=['Header',''])
        # df_pressure.columns = columns_pressure

        df_pressure = df_mz['TP99']['Ion Current']
  
    
        # Machine Time Data
        df_machine_time = pd.read_csv(Filename_MS, skiprows=7, sep='\t', decimal=',').iloc[:, 1]
        machine_time_header = pd.read_csv(Filename_MS, skiprows=5, nrows=0, sep='\t\t\t', engine='python').columns.to_list()[:0]
        columns_machine_time = pd.MultiIndex.from_product([machine_time_header, ['time[min]']], names=['Header',''])
        df_machine_time.columns = columns_machine_time
    
        mz_data = {'m/z' : df_mz, 'pressure' : df_pressure, 'time' : df_machine_time}
        # print(mz_data['pressure'])
        # print((df_mz['44']['Ion Current'])/mz_data['pressure'])
        self.mz_data = mz_data


        # default_Temp = r'C:\Users\bjorngso\OneDrive - Universitetet i Oslo\01 Results\Testing\Temperature\SW3000FC_5-26-2021\T3000FC_26_5_2021-13-2-18-869.csv'
        # Filename_Temp = kwargs.get('file_Temp', default_Temp)
        
        # header_temp = pd.read_csv(Filename_Temp, skiprows=9,nrows=0,sep=';',engine='python')
        # df_temp =pd.read_csv(Filename_Temp, skiprows=10, sep=';', decimal=',',engine='python', header=None)

        # header_temp = header_temp.columns.to_list()
        # header_new = [f[:-1] for f in header_temp]
        # df_temp.columns = header_new
        
        # time = pd.to_datetime(df_temp['Max Time'].values, format="%d.%m.%Y %H:%M:%S")
        # time_d = (time - time[0])/np.timedelta64(1, 'm')
        # self.time_d = time_d
        # df_temp['Max Time'] = time_d
        # self.df_temp=df_temp
        return 

    
    

#############################################################################################
#                                                                                           #
#                                          Misc                                             #
#                                                                                           #
#                                                                                           #
#############################################################################################

    def give_values(self, **kwargs):
        Temp_rate,subtracted_CO2 = self.CO2_fitted_XY
        Temp_rate,subtracted_CH4 = self.CH4_fitted_XY
        CO2_He = np.array(self.CO2_He)
        He_CH4 = np.array(self.He_CH4)
        T_norm = np.array(self.T_norm)
        return CO2_He,He_CH4,T_norm, subtracted_CO2,subtracted_CH4,Temp_rate






#############################################################################################
#                                                                                           #
#                                   Normal MTM analysis                                     #
#                                   -------------------                                     #
#                                                                                           #
#############################################################################################





    def plotting_MS(self,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12, fig, ax, peak_start, peak_end,plot_start_end):

        ax.clear()
        mz_data = self.mz_data
        df_mz = df = (mz_data['m/z'])
        T = (mz_data['time']/60).to_numpy()
        P = (mz_data['pressure']).to_numpy()

        H= (df_mz['2','Ion Current']/P).to_numpy()
        MeOH_29= (df_mz['29','Ion Current']/P).to_numpy()
        MeOH_31= (df_mz['31','Ion Current']/P).to_numpy()
        CO2 = (df_mz['44','Ion Current']/P).to_numpy()
        DME = (df_mz['46','Ion Current']/P).to_numpy()
        H2O = (df_mz['18','Ion Current']/P).to_numpy()
        CO = (df_mz['28','Ion Current']/P).to_numpy()
        CH4_16 = (df_mz['16','Ion Current']/P).to_numpy()
        CH4_15 = (df_mz['15','Ion Current']/P).to_numpy()
        He = (df_mz['4','Ion Current']/P).to_numpy()
        Ne = (df_mz['20','Ion Current']/P).to_numpy()
        
        if w1 is True:
            ax.plot(T, MeOH_29, label='MeOH (29)')
        if w2 is True:
            ax.plot(T, MeOH_31, label='MeOH (31)')
        if w3 is True:
            ax.plot(T, CO2, label='CO2')
        if w4 is True:
            ax.plot(T, DME, label='DME')
        if w5 is True:
            ax.plot(T, CO, label='CO')
        if w6 is True:
            ax.plot(T, H2O, label='H2O')
        if w7 is True:
            ax.plot(T, CH4_15, label='CH4_15')
        if w8 is True:
            ax.plot(T, CH4_16, label='CH4_16')
        if w9 is True:
            ax.plot(T, He, label='He')
        if w11 is True:
            ax.plot(T, Ne, label='Ne')
        if w12 is True:
            ax.plot(T, H, label='H')
            

        ax.legend()

        if w10 == 'Logarythmic':
            ax.set_yscale('log')
        else: pass

        ax.grid(True)
        ax.set_title('MS data')
        ax.set_xlabel('Extraction time [$min$]')
        ax.set_ylabel('Ion Current [$m/z$]')
        
        
        if 'peak_start'  in plot_start_end.keys():
            del plot_start_end['peak_start']
            
        elif 'peak_end'  in plot_start_end.keys():
            del plot_start_end['peak_end']
            
        plot_start_end.update({'peak_start': peak_start})
        plot_start_end.update({'peak_end': peak_end})
        self.plot_start_end = plot_start_end
    
        return


    def plot_MS(self,**kwargs):
        fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=False)
        plot_start_end = {}   
        w = widgets.interactive(self.plotting_MS,
                w1 = widgets.Checkbox(value= True, description='MeOH_29',disabled=False,indent=False),
                w2 = widgets.Checkbox(value= False, description='MeOH_31',disabled=False,indent=False),
                w3 = widgets.Checkbox(value= False, description='CO2',disabled=False,indent=False),
                w4 = widgets.Checkbox(value= False, description='DME',disabled=False,indent=False),
                w5 = widgets.Checkbox(value= False, description='CO',disabled=False,indent=False),
                w6 = widgets.Checkbox(value= False, description='H2O',disabled=False,indent=False),
                w7 = widgets.Checkbox(value= False, description='CH4_15',disabled=False,indent=False),
                w8 = widgets.Checkbox(value= False, description='CH4_16',disabled=False,indent=False),
                w9 = widgets.Checkbox(value= False, description='He',disabled=False,indent=False),
                w10 = widgets.Dropdown(options=['Logarythmic','Linear'], description='y-axis format'),
                w11 = widgets.Checkbox(value= False, description='Ne',disabled=False,indent=False),
                w12 = widgets.Checkbox(value= False, description='H',disabled=False,indent=False),
                fig=widgets.fixed(fig), ax=widgets.fixed(ax),
                peak_start = widgets.FloatText(value='1200', placeholder='x-value', description='Start x value', disabled=False,layout=widgets.Layout(width='20%')),
                peak_end = widgets.FloatText(value='1350', placeholder='x-value', description='End x value', disabled=False,layout=widgets.Layout(width='20%')),
                plot_start_end=widgets.fixed(plot_start_end),
                )  
        return display(w)


#############################################################################################
#                                                                                           #
#                                     define the Area                                       #
#                                                                                           #
#                                                                                           #
#############################################################################################

    def Area(self,start,stop, T, MS):

        id_min = np.where((T>start) & (T < (start+0.1)))[0][0]
        id_max = np.where((T>stop) & (T < stop+0.1))[0][0]

        T_select = np.array([T[id_max], T[id_min]])
        MS_select = np.array([MS[id_max], MS[id_min]])

        a = np.diff(MS_select)/np.diff(T_select) # a
    #    y= MS[id_max]-(a*T[id_max])
        x= MS[id_min]-(a*T[id_min])
        
        area_MS=sp.trapz(self.f_cor_MS(T[id_min:id_max], MS[id_min:id_max], a, x), T[id_min:id_max],)
        baseline = self.f_baseline(T, a, x)
        
        return baseline[id_min:id_max], T[id_min:id_max], format(float(area_MS),'.5e'), MS[id_min:id_max] 

    def f_baseline(self,T, a, x):
        y_line = a*T+x
        return y_line

    def f_MS(self,MS):
        return MS
    
    def f_cor_MS(self,T, MS, a, x):
        y_cor_MS= self.f_MS(MS)-self.f_baseline(T, a, x)
        return y_cor_MS





    def plotting_Area(self,start, stop, mz, fig, ax, Click,Data_for_prossecing):
    
        mz_dict = {
            'CO2' : '44',
            'DME' : '46',
            'MeOH_29' : '29',
            'MeOH_31' : '31',
            'CO' : '28',
            'H2O' : '18',
            'CH4_16' : '16',
            'CH4_15' : '15'
        }
        

        mz_data = self.mz_data
        df_mz = (mz_data['m/z'])
        T = (mz_data['time']/60).to_numpy()
        P = (mz_data['pressure']).to_numpy()
        mz_ion = (df_mz[mz_dict[mz]]['Ion Current']/P).to_numpy()
        
        plot_start_end = self.plot_start_end

        ax.clear()

        

       
        ion_baseline, T_ion, ion_area, ion_signal = self.Area(start, stop, T, mz_ion)
        ax.plot(T_ion, ion_baseline, label='baseline', color='black')
        ax.plot(T, mz_ion, label=mz,color='red')
        ax.set_ylim(0,(np.max(mz_ion[:-1])))

        inx_min = np.where((T>plot_start_end['peak_start']) & (T<plot_start_end['peak_end']))[0][0]
        inx_max = np.where((T>plot_start_end['peak_start']) & (T<plot_start_end['peak_end']))[0][-1]
        
        max_peak = np.max(mz_ion[inx_min:inx_max])
        idx = np.where(mz_ion == max_peak)[0]

        ax.set_ylim(0,(np.max(mz_ion[idx])+abs(np.max(np.diff(mz_ion[inx_min:inx_max])))))
      

                
        ax.grid(True)
        ax.set_title('MS data integrated')
        ax.set_xlabel('Time [$min$]')
        ax.set_ylabel('Ion Current [$m/z$]')
        ax.fill_between(T_ion,ion_baseline, ion_signal, color='lightyellow', label='Integrated area')
        ax.set_xlim(plot_start_end['peak_start'], plot_start_end['peak_end'])
        ax.legend() 
        
        
        ax.text(T_ion.mean(), 3.0e-6 , '{}'.format(ion_area,'.5e'), bbox={'facecolor': 'lightgreen', 'alpha': 1, 'pad': 5})
        print('The integrated area is:', ion_area)

        
        
        if Click is True:
            if mz in Data_for_prossecing: 
                del Data_for_prossecing[mz]
            Data_for_prossecing.update({mz: float(ion_area)})

        return Data_for_prossecing
         



    def plot_MS_Area(self,**kwargs):
        fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=False)
        Data_for_prossecing = {}
        
            
        w = widgets.interactive(self.plotting_Area, start=(self.plot_start_end['peak_start'], (self.plot_start_end['peak_end']-1), 0.01), stop=((self.plot_start_end['peak_start']+1),  self.plot_start_end['peak_end'], 0.01), mz=widgets.Dropdown(
            options=['CO2', 'DME','MeOH_29','MeOH_31','CO','H2O','CH4_16','CH4_15'],
            value='CO2',
            description='Mass:',
            disabled=False,
            ),
            fig=widgets.fixed(fig), ax=widgets.fixed(ax),
            Click = widgets.Checkbox(value=False, description='Save Area and add to dictionary',disabled=False,button_style='', tooltip='Description', icon='check'),
            Data_for_prossecing = widgets.fixed(Data_for_prossecing) 
            )
        self.Data_for_prossecing = Data_for_prossecing
        return display(w)




    def Area_to_csv(self,**kwargs):
        filename = self.Filename_MS
        Data_for_prossecing = self.Data_for_prossecing
        filename_truncated = filename[82:-4]
        print('{}_Areas.csv'.format(filename_truncated))
        Data = pd.DataFrame.from_dict(Data_for_prossecing,orient='index')

        if os.path.exists('{}_Areas.csv'.format(filename_truncated)):
            os.remove('{}_Areas.csv'.format(filename_truncated))
            
        Data.to_csv('{}_Areas.csv'.format(filename_truncated))

    def where_am_i(self):
        print(os.getcwd())
        



#############################################################################################
#                                                                                           #
#                         Calculating test data from raw data                               #
#                                                                                           #
#                                                                                           #
#############################################################################################


    def moles(self,Calib, Area, Flow, Pressure, R, Temperature):
        if Area is False:
            Area = 0
        def volume(Calib, Area, Flow):
            def percentage(Calib,Area):
                percent = Calib*(100/1000000) #    % min
                return percent
            vol = (percentage(Calib, Area)/100)*Flow #    mL
            return vol
        mol = (float(Pressure)*(volume(float(Calib), float(Area), float(Flow))*10e-3))/(float(R)*float(Temperature)) #mol

        return mol

    
    def calib_210908(self, x, type):
        #from calibration of the MS 20.07.2021

        if type == "MeOH_29":
            y = x/(6.904*10**(-8))
        elif type == "MeOH_31":
            y = (x+(5*10**(-6)))/(2*10**(-8))
        elif type == "DME":
            y = x/(4.337*10**(-8))
        elif type == "CO2":
            y = x/(9.955*10**(-8))
        elif type == "CO":
            y = x/(4*10**(-7))
        else:
            pass
        return y # mol ppm

    def new_MTM_210422(self, x, type):
        #from calibration of the MS 20.07.2021

        if type == "MeOH_29":
            y = 4*10**(6*x)
        elif type == "MeOH_31":
            y = -2*10**(10*x**2)+2*10**(7*x)
        elif type == "DME":
            y = 10**(7*x)
        elif type == "CO2":
            y = 3*10**(6*x)
        elif type == "CO":
            y = 916512*x
        else:
            pass
        return y # mol ppm

    def new_MTM_220515(self, x, type):
        Calibration = {'CO' : 1, 
                        'CO2' : 2.917*10**(-8),
                        'DME' : 1.174*10**(-8), 
                        'MeOH_29' : 1,  
                        'MeOH_31' : 9.434*10**(-9)}
        consentration = x/Calibration[type] #ppm min
        return consentration
    

    def new_HP_MTM_2022(self, x, type):
        Calibration = {'CO' : 1,
                       'CO2' : 3.686*10**(-7),
                       'DME' : 9.804*10**(-8),
                       'MeOH_29' : 2.707*10**(-7),
                       'MeOH_31' : 9.491*10**(-8)}
        consentration = x/Calibration[type] #ppm min
        return consentration

    def calib_old(self,x,type):
        Calibration = {'CO' : self.CO_calib_factor, 
                        'CO2' : self.CO2_44_calib_factor,
                        'DME' : self.DME_46_calib_factor, 
                        'MeOH_29' : self.MeOH_29_calib_factor, 
                        'MeOH_31' : self.MeOH_31_calib_factor}
        consentration = x/Calibration[type] #ppm min
        return consentration


    def Cu_calc(self,SiAl,CuAl,SAPO):
        Al = 26.98                                            #g/mol
        P = 30.97                                             #g/mol
        Si = 28.09                                            #g/mol
        O = 16                                                #g/mol
        Cu = 63.546                                           #g/mol
        PAlav = 28.975                                        #g/mol

        if SAPO == True:
            AlP_Si = SiAl
            Cu_Si = CuAl
            Mm_empty = (AlP_Si*PAlav)+(1*Si)+(2*(Cu_Si+1)*O)  #g/mol
            Mm_loaded = Mm_empty+Cu_Si*Cu                     #g/mol with Cu
            wt = 100*(Cu_Si*Cu)/Mm_loaded                     #wt%
            micromol_g = 1000000*(wt/(100*Cu))                #micro mol/g
            
        else:
            Si_Al = SiAl
            Cu_Al = CuAl
            Mm_empty = (Si_Al*Si)+(1+Al)+(2*(Si_Al+1)*O)     #g/mol
            Mm_loaded = Mm_empty+Cu_Al*Cu                    #g/mol with Cu
            wt = 100*(Cu_Al*Cu)/Mm_loaded                    #wt%
            micromol_g = 1000000*(wt/(100*Cu))               #micro mol/g
            
        return micromol_g






    def manual_inputting(self,**kwargs):
        file_name_ID = kwargs.get('file_name_ID','')
       
        Mass_cat = kwargs.get('M_c',0.0)
        Water_percentage = kwargs.get('W_p',0.0)
        Dry_weight = kwargs.get('D_w',0.0)
        SiAl = kwargs.get('SiAl',0.0)
        CuAl = kwargs.get('CuAl',0.0)
        Cu_content = kwargs.get('Cu_content',0.0)
        SAPO = kwargs.get('SAPO',True)
        all_input = kwargs.get('all_input',False)
        csv_file = kwargs.get('csv_file','')

        Temperature = 298.15 #K
        R = 0.082057 #L atm K-1 mol-1
        Cu = self.Cu_calc(float(SiAl),float(CuAl),SAPO)
        Flow = kwargs.get('F',0)   
        Pressure = kwargs.get('P',0)
        calibration_using = kwargs.get('calibration_using', 'old')

        if all_input is True:

            CSV_file = r'C:\Users\bjorngso\OneDrive - Universitetet i Oslo\03 Codes\Jupter-Data\Data analysis\ACTIVE_SCRIPTS\{}'.format(csv_file)
            Area = self.Read_MZ_Area_from_csv(CSV_file)

            for key in Area:
                if Area[key] < 0:
                    Area[key] = 0
                    print(f'-   The area for {key} is negative, and will be sett to {0}.')
            if 'CO' not in Area.keys():
                Area['CO'] = 0
                print(f'-   CO is not in dict, and will be added.')
            if 'DME' not in Area.keys():
                Area['DME'] = 0
                print(f'-   DME is not in dict, and will be added.')
        
    
            Calibration = {'CO' : self.CO_calib_factor, 
                        'CO2' : self.CO2_44_calib_factor,
                        'DME' : self.DME_46_calib_factor, 
                        'MeOH_29' : self.MeOH_29_calib_factor, 
                        'MeOH_31' : self.MeOH_31_calib_factor}

            Yield = {}
            
            for key in Calibration:
                
                if calibration_using == 'old':
                    moles = format(self.moles(self.calib_old(float(Area[key]),key),float(Area[key]),Flow,Pressure,R,Temperature), '.2e')
                elif calibration_using =='210908':
                    moles = format(self.moles(self.calib_210908(float(Area[key]),key),float(Area[key]),Flow,Pressure,R,Temperature), '.2e')
                elif calibration_using =='210422':
                    moles = format(self.moles(self.new_MTM_210422(float(Area[key]),key),float(Area[key]),Flow,Pressure,R,Temperature), '.2e')
                elif calibration_using =='220515':
                    moles = format(self.moles(self.new_MTM_220515(float(Area[key]),key),float(Area[key]),Flow,Pressure,R,Temperature), '.2e')
                elif calibration_using == 'HP':
                    moles = format(self.moles(self.new_HP_MTM_2022(float(Area[key]),key),float(Area[key]),Flow,Pressure,R,Temperature), '.2e')
                Yields = np.round(((float(moles)*10**5)/float(Dry_weight)), decimals=1)
                Yield["%s" %key] = Yields

            
            Yield_MeOH_29 = np.round(Yield['MeOH_29']+(2*Yield['DME']))
            Yield_MeOH_31 = np.round(Yield['MeOH_31']+(2*Yield['DME']))
            
            Selectivity = np.round(Yield_MeOH_31/(Yield['CO2']+Yield_MeOH_31+Yield['CO']),decimals=2)

            if Cu == 0 or Cu == None:
                if Cu_content == 0 or Cu_content == None:
                    Productivity = 0
                    Productivity31 = 0
                else: 
                    Cu = Cu_content
                    Productivity = Yield_MeOH_29/Cu
                    Productivity31 = Yield_MeOH_31/Cu
            else:
                Productivity = Yield_MeOH_29/Cu
                Productivity31 = Yield_MeOH_31/Cu
                
            D = [{'micro mol per gram zeolite':'',
                    'CO yield':Yield['CO'],
                    'MeOH yield':Yield['MeOH_29'],
                    'MeOH(31) yield':Yield['MeOH_31'],
                    'DME yield':Yield['DME'],
                    'CO2 yield':Yield['CO2'],
                    'Yield':Yield_MeOH_29, 
                    'Yield(31)':Yield_MeOH_31, 
                    'Selectivity':Selectivity, 
                    'Productivity':Productivity,
                    'Productivity(31)':Productivity31}]
            Df= pd.DataFrame.from_dict(D)

            print('')
            display(Df)


            sample_ID = file_name_ID
            if sample_ID == None or sample_ID == '':
                print('')
                print('')
                print('')
                print('')
                print('')
                print('                  ----------   No excel file name wass given, and a file will not be created    ----------')
            else:
                if os.path.exists(r"C:\Users\bjorngso\OneDrive - Universitetet i Oslo\01 Results\Testing\M2M\{}_results.xlsx".format(sample_ID)):
                    os.remove(r"C:\Users\bjorngso\OneDrive - Universitetet i Oslo\01 Results\Testing\M2M\{}_results.xlsx".format(sample_ID))

                Df.to_excel(r"C:\Users\bjorngso\OneDrive - Universitetet i Oslo\01 Results\Testing\M2M\{}_results.xlsx".format(sample_ID))
                self.DF = Df
        return






    def filen_name_excel(self,**kwargs):
        file_name_ID = kwargs.get('file_name_ID',[])
        file_name_lst = kwargs.get('file_name_lst',[])
        file_name_lst.append(file_name_ID)
        return

    def Read_MZ_Area_from_csv(self,Filename):
        df_Area = pd.read_csv(Filename, skiprows=1,header=None, names =['Mass','Area'], index_col=[0], usecols = ['Mass','Area'])
        dict_Area = pd.DataFrame.to_dict(df_Area)['Area']
        return dict_Area





    def testing_resutls(self,**kwargs):

        w = widgets.interactive(self.manual_inputting,
                                csv_file = widgets.Text(value=r'TEST_CSV.csv', placeholder='Sample are values.', description=r'CSV-file:', disabled=False,layout=widgets.Layout(width='90%')),
                                calibration_using = widgets.Dropdown(options=[('Autumn 2020 (old MTM)','old'),('Summer 2021','210908'),('spring 2021','210422'),('spring 2022 (new MTM)','220515'),('HP_MTM','HP')],description='Calibration',disabled=False,layout=widgets.Layout(width='20%')),
                                F = widgets.FloatText(value='16.5', placeholder='Flow, during extraction', description=r'F (mL/min):', disabled=False,layout=widgets.Layout(width='20%')),
                                P = widgets.FloatText(value='1.0', placeholder='Pressure, during extraction', description=r'P (atm):', disabled=False,layout=widgets.Layout(width='20%')),
                                M_c = widgets.FloatText(value='0.1', placeholder='Mass catalyst', description=r'm_cat (g):', disabled=False,layout=widgets.Layout(width='20%')),
                                W_p = widgets.FloatText(value='10.0', placeholder='Water percentage in seolite', description='H2O %:', disabled=False,layout=widgets.Layout(width='20%')),
                                D_w = widgets.FloatText(value='0.086', placeholder='Weight cat, after TGA', description=r'm_TGA (g):', disabled=False,layout=widgets.Layout(width='20%')),
                                SiAl = widgets.FloatText(value='0', placeholder='Si/Al:', description=r'Si/Al::', disabled=False,layout=widgets.Layout(width='20%')),
                                CuAl = widgets.FloatText(value='0', placeholder='Cu/Al Ratio', description=r'Cu/Al:', disabled=False,layout=widgets.Layout(width='20%')),
                                Cu_content = widgets.FloatText(value='0', placeholder='Cu_content', description=r'Cu_content:', disabled=False,layout=widgets.Layout(width='20%')),
                                SAPO = widgets.Checkbox(value=False,description='Is the cat SAPO?',disabled=False,indent=False),
                                file_name_ID = widgets.Text(value=r'', placeholder='data_name_experiman-type (2021-04-12_0-11CuMOR6-5_3)', description='Result file name:', disabled=False,layout=widgets.Layout(width='90%')),
                                all_input = widgets.Checkbox(value=False,description='Is all input given?',disabled=False,indent=False),
                                )

        return display(w)


















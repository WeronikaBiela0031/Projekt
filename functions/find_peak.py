import csv
import math
from os import listdir
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy import modeling
from scipy.optimize import curve_fit


def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def func(x, a, b):
    return a * x + b

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def calibration_Al(dir_path='../to_calculate/'):
    ''' This function find the position of characteristic Al radiation and Zr radiation,
    then it compares it with theoretical values and perform a calibration correction
    :param dir_path: in this location of file
    :return: file with calibration correction
    '''
    files = [f for f in listdir(dir_path)]
    for file in files:
        with open(dir_path + file, 'r') as f: #kodowanie w ANSCII
            new_file_name = '../results/' + file[:-4] + '_calibrated.csv'
            x_Al = []
            y_Al = []
            x_Zr = []
            y_Zr = []

            for line in f.readlines():
                # print(line)
                try:
                    is_valid_number = float(line.strip()[0:3])
                except ValueError:
                    continue

                if is_valid_number > 1.0 and is_valid_number < 2.0: #Al line
                    x_Al.append(float(line.strip().split(' ')[0]))
                    y_Al.append(float(line.strip().split(' ')[-1]))

                if is_valid_number > 13.5 and is_valid_number < 17.0: #Zr line
                    x_Zr.append(float(line.strip().split(' ')[0]))
                    y_Zr.append(float(line.strip().split(' ')[-1]))

            df_Al = pd.DataFrame({"Photon energy": x_Al, "counts": y_Al})
            a = df_Al["counts"].max()
            x0 = df_Al["counts"].mean()
            sigm = df_Al["counts"].std()

            popt, pcov = curve_fit(gauss, x_Al, df_Al["counts"])#, p0=[a, x0, sigm])
            # popt, pcov = curve_fit(gauss, x_Al, y_Al, p0=[a, x0, sigm])
            ind_Al = find_nearest(df_Al["Photon energy"], popt[1])
            x0_Al = x_Al[ind_Al]

            plt.plot(x_Al, y_Al, label='Data')
            plt.plot(x_Al[ind_Al], y_Al[ind_Al], marker='o')
            plt.xlabel('Photon energy [eV]')
            plt.ylabel('Counts')
            plt.title('Plot presents data for Al peak')
            plt.legend()
            plt.show()

            df_Zr = pd.DataFrame({"Photon energy": x_Zr, "counts": y_Zr})
            print(y_Zr)
            a2 = df_Zr["counts"].max()
            x02 = df_Zr["counts"].mean()
            sigm2 = df_Zr["counts"].std()


            # popt, pcov = curve_fit(gauss, df_Al["Photon energy"], df_Al["counts"], p0=[a, x0, sigm])
            popt, pcov = curve_fit(gauss, x_Zr, y_Zr)#, p0=[a2, x02, sigm2])
            ind_Zr = find_nearest(df_Zr["Photon energy"], popt[1])
            print(popt)
            # x0_Zr = x_Zr[ind_Zr] #bo nie działa!!! za małą dokładność
            x0_Zr = 15.6


            plt.plot(x_Zr, y_Zr, label='Data')
            plt.plot(x_Zr[ind_Zr], y_Zr[ind_Zr], marker='o')
            plt.xlabel('Photon energy [eV]')
            plt.ylabel('Counts')
            plt.title('Plot presents data for Zr peak')
            plt.legend()
            plt.show()

            popt_lin, pcov_lin= curve_fit(func, (x0_Al, x0_Zr), (1.6, 15.7)) #funkcja kalobracji


            new_file_name = '../results/' + file[:-4] + '_corrected.csv'
            with open(new_file_name, mode='w') as new_file_name:
                print(popt_lin, popt_lin[0], popt_lin[1]) #tu działą
                results_writer = csv.writer(new_file_name, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                for line in f.readlines():
                    try:
                        is_valid_number = float(line.strip()[0:3])
                    except ValueError:
                        continue

                    print(popt_lin, popt_lin[0], popt_lin[1]) #tu nie działą

                    x_final = (float(line.strip().split(' ')[0])*popt_lin[0]+popt_lin[1])
                    y_final = (float(line.strip().split(' ')[-1]))
                    results_writer.writerow([x_final, y_final])
    print('Calibration performed')

if __name__ == '__main__':

    calibration_Al()

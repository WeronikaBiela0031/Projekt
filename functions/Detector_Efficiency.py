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


def round_half_up(n):  # w przypadku naszych obliczeń chcemy zaokrąglić eV do całości
    return math.floor(n + 0.5)


def detector_eff(dir_path='../to_calculate/'):
    ''' This program is providing the detector efficiency correction caused by absorption of 25 um layer of beryllium
    :param a path to catalog in which we have saved the data for detector efficiency correction
    :return: function does not return anything
    '''
    files = [f for f in listdir(dir_path)]
    dict_file_Be = {}
    with open('25uBe.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            dict_file_Be[float(row[0])] = float(row[1])

    for file in files:
        with open(dir_path + file, 'r') as f:
            new_file_name = '../results/' + file[:-4] + '_corrected.csv'
            x = []
            y = []
            y_new = []
            with open(new_file_name, mode='w') as new_file_name:
                results_writer = csv.writer(new_file_name, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                y_max = 0  # potrzebny limit dla liczby zliczeń na wykresie
                for line in f.readlines():
                    try:
                        is_valid_number = float(line.strip()[0:3])
                    except ValueError:
                        continue
                    else:
                        if is_valid_number < 0.5 or is_valid_number > 19:
                            continue
                    photon_energy = round_half_up(float(line.strip().split(' ')[0]) * 1000)
                    coefficient = dict_file_Be.get(photon_energy)
                    if coefficient == 0:
                        continue
                    counts_raw = float(line.strip().split(' ')[-1])
                    counts = counts_raw / coefficient
                    results_writer.writerow([photon_energy, counts])
                    if photon_energy > 2500 and y_max < counts:
                        y_max = counts

                    x.append(photon_energy)
                    y.append(counts_raw)
                    y_new.append(counts)

                plt.plot(x, y, label='Row data')
                plt.xlabel('Photon energy [eV]')
                plt.ylabel('Counts')
                plt.title('Plot presents data before absorption correction')
                plt.legend()
                plt.show()

                plt.xlim(2500, 20000)
                plt.ylim(-500, y_max)
                plt.plot(x, y_new, label='Corrected data')
                plt.xlabel('Photon energy [eV]')
                plt.ylabel('Counts')
                plt.title('Plot presents data after absorption correction')
                plt.legend()
                plt.show()

    print("Detector Efficiency Applied")


def detector_eff_aluminium(dir_path='../to_calculate/'):
    ''' This program is providing the detector efficiency correction caused by absorption
    of 25 um layer of beryllium and 70 um layer of Aluminium (if it is necessary)
    :param a path to catalog in which we have saved the data for detector efficiency correction
    :return: function does not return anything
    '''
    files = [f for f in listdir(dir_path)]
    dict_file_Be = {}
    dict_file_Al = {}

    with open('70uAl.csv') as csv_file_Al:
        csv_reader_Al = csv.reader(csv_file_Al, delimiter=',')
        for row in csv_reader_Al:
            dict_file_Al[float(row[0])] = float(row[1])
    with open('25uBe.csv') as csv_file_Be:
        csv_reader_Be = csv.reader(csv_file_Be, delimiter=',')
        for row in csv_reader_Be:
            dict_file_Be[float(row[0])] = float(row[1])

    for file in files:
        if file.count('Al70') == 0:
            D = dict_file_Al.copy()
            dict_file_Al.clear()
            for key in D.keys():
                dict_file_Al.setdefault(key, 1)

        with open(dir_path + file, 'r') as f:
            new_file_name = '../results/' + file[:-4] + '_Al_corrected.csv'
            with open(new_file_name, mode='w') as new_file_name:
                results_writer = csv.writer(new_file_name, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                for line in f.readlines():
                    try:
                        is_valid_number = float(line.strip()[0:3])
                    except ValueError:
                        continue
                    else:
                        if is_valid_number < 0.9 or is_valid_number > 19:
                            continue
                    photon_energy = round_half_up(float(line.strip().split(',')[0]) * 1000)

                    coefficient = dict_file_Be.get(photon_energy) * dict_file_Al.get(photon_energy)
                    if coefficient == 0:
                        continue
                    try:
                        counts_raw = float(line.strip().split(',')[-1])
                    except ValueError:
                        counts_raw = 0
                    counts = counts_raw / coefficient
                    results_writer.writerow([photon_energy, counts])
    print("Correction for both beryllium and aluminium window has been applied")


def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def calibration_Al(dir_path='../to_calculate/'):
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
            x0_Zr = x_Zr[ind_Zr]

            plt.plot(x_Zr, y_Zr, label='Data')
            plt.plot(x_Zr[ind_Zr], y_Zr[ind_Zr], marker='o')
            plt.xlabel('Photon energy [eV]')
            plt.ylabel('Counts')
            plt.title('Plot presents data for Zr peak')
            plt.legend()
            plt.show()
    print('Calibration performed')

if __name__ == '__main__':
   # detector_eff()
    #detector_eff_aluminium()
    calibration_Al()

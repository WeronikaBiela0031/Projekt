import csv
import math
from os import listdir
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy import optimize


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


def _1gaussian(x, amp1, cen1, sigma1):
    return amp1 * (1 / (sigma1 * (np.sqrt(2 * np.pi)))) * (
        np.exp((-1.0 / 2.0) * (((x - cen1) / sigma1) ** 2)))

def calibration_Al(dir_path='../to_calculate/'):
    files = [f for f in listdir(dir_path)]
    for file in files:
        with open(dir_path + file, 'r') as f: #kodowanie w ANSCII
            new_file_name = '../results/' + file[:-4] + '_calibrated.csv'
            x = []
            y = []
            y_new = []
            # with open(new_file_name, mode='w') as new_file_name:
            #     results_writer = csv.writer(new_file_name, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for line in f.readlines():
                # print(line)
                try:
                    is_valid_number = float(line.strip()[0:3])
                except ValueError:
                    continue

                if is_valid_number > 1.0 and is_valid_number < 2.0: #Al line

                    x.append(line.strip().split(' ')[0])
                    y.append(line.strip().split(' ')[-1])
            df_Al = pd.DataFrame({"Photon energy": x, "counts": y})
            import pdb; pdb.set_trace()
            popt_gauss, pcov_gauss = scipy.optimize.curve_fit(_1gaussian, x, y,
                                                              p0=[amp1, cen1, sigma1])
            perr_gauss = np.sqrt(np.diag(pcov_gauss))
            peaks, properties = find_peaks(y, width= 20)
            plt.plot(y)
            plt.plot(peaks, y[peaks], "counts")
            plt.show()
    print('Calibration performed')

if __name__ == '__main__':
   # detector_eff()
    #detector_eff_aluminium()
    calibration_Al()

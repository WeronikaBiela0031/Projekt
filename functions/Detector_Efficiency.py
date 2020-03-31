import csv
from os import listdir


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
            with open(new_file_name, mode='w') as new_file_name:
                results_writer = csv.writer(new_file_name, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                for line in f.readlines():
                    try:
                        is_valid_number = float(line.strip()[0:3])
                    except ValueError:
                        continue
                    else:
                        if is_valid_number < 0.5 or is_valid_number > 19:
                            continue
                    # TODO fix rounding numbers
                    photon_energy = round(float(line.strip().split(' ')[0]) * 1000)
                    coefficient = dict_file_Be.get(photon_energy)
                    if coefficient == 0:
                        continue
                    counts_raw = float(line.strip().split(' ')[-1])
                    counts = counts_raw / coefficient
                    results_writer.writerow([photon_energy, counts])

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
            new_file_name = '../results/' + file[:-4] + '_corrected.csv'
            with open(new_file_name, mode='w') as new_file_name:
                results_writer = csv.writer(new_file_name, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                for line in f.readlines():
                    try:
                        is_valid_number = float(line.strip()[0:3])
                    except ValueError:
                        continue
                    else:
                        if is_valid_number < 0.5 or is_valid_number > 19:
                            continue
                    # TODO fix rounding numbers
                    photon_energy = round(float(line.strip().split(' ')[0]) * 1000)
                    coefficient = dict_file_Be.get(photon_energy)*dict_file_Al.get(photon_energy)
                    if coefficient == 0:
                        continue
                    counts_raw = float(line.strip().split(' ')[-1])
                    counts = counts_raw / coefficient
                    results_writer.writerow([photon_energy, counts])

    print("Correction for aluminium window has been applied")


if __name__ == '__main__':
    detector_eff()
   # detector_eff_aluminium()

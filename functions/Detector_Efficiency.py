import csv
from os import listdir


def detector_eff(dir_path='../to_calculate/'):
    ''' This program is providing the detector efficiency correction caused by absorption of 25 um layer of beryllium
    :param a path to catalog in which we have saved the data for detector efficiency correction
    :return:
    '''
    files = [f for f in listdir(dir_path)]
    dict_file = {}
    with open('25uBe.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            dict_file[float(row[0])] = float(row[1])

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
                    coefficient = dict_file.get(photon_energy)
                    if coefficient == 0:
                        continue
                    counts_raw = float(line.strip().split(' ')[-1])
                    counts = counts_raw / coefficient
                    results_writer.writerow([photon_energy, counts])

    print("Detector Efficiency Applied")


if __name__ == '__main__':
    detector_eff()

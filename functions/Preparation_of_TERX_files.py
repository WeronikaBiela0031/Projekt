# TODO 1 rzutowanie wszystkich kolumn,
#     2 zczytanie kalibracji wynikającej z pliku par
#     3 Spectrum_Calibration
import json
from os import listdir


def projection_of_TERX_file(dir_path='../data/terx/', start=[0],end=[-1]):  # time_interval domyslnie od początku- index 0 do końca index -1
    files_list = [f for f in listdir(dir_path) if f[-3:] != 'par']
    files_par = [f for f in listdir(dir_path) if f[-3:] == 'par']
    files_dict = dict(zip(files_list, files_par))
    counts_vs_photon_energy = {}
    for f in files_list:
        with open(dir_path + f, 'r') as my_file, open(dir_path + files_dict[f]) as my_par_file:
            par_line = my_par_file.readlines()
            E0, E_bin, t_bin, t_start = [float(x) for x in par_line[0].strip("\n").split("\t")]
            E0 = E0 * 1000  # eV

            lines = my_file.readlines()
            for line in lines:
                counts_vs_photon_energy[E0] = sum([int(x) for x in line.strip('\n').split('\t')])
                E0 += float(E_bin)
                # import pdb; pdb.set_trace()
            with open('testowy_rezultat.json', 'w') as t:
                t.write(json.dumps(counts_vs_photon_energy))

    print('Projection done')

# TODO

if __name__ == '__main__':
    projection_of_TERX_file()

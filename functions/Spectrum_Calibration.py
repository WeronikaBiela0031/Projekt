from os import listdir

import pandas as pd
from scipy.signal import find_peaks


def calibration(dir_path='../data/'):
    ''' This function takes a TERX file and a corresponding SPX file (save in a ANSI encoding) and performs the calibration
    :param dir_path: in this location there needs to be two catalogs named: terx and spx, consisting respectively proper files
    :return:file with a parameters with a correct calibration for TERX file
    '''
    files_list = [f for f in listdir(dir_path + "terx/") if f[-3:] != 'par']
    files_par = [f for f in listdir(dir_path + "terx/") if f[-3:] == 'par']
    files_spx = [f for f in listdir(dir_path + "spx/")]
    files_dict = dict(zip(files_list, files_par))
    counts_vs_photon_energy = {}
    for f in files_list:
        with open(dir_path + "terx/" + f, 'r') as my_file, open(dir_path + 'terx/' + files_dict[f]) as my_par_file:
            for spx in files_spx:
                if f[:8] == spx[:8]:
                    with open(dir_path + 'spx/' + spx) as my_spx:
                        par_line = my_par_file.readlines()
                        E0_first, E_bin, t_bin, t_start = [float(x) for x in par_line[0].strip("\n").split("\t")]
                        E0 = E0_first * 1000  # eV

                        lines = my_file.readlines()
                        for line in lines:
                            counts_vs_photon_energy[E0] = sum([int(x) for x in line.strip('\n').split('\t')])
                            E0 += float(E_bin)
                        df_terx = pd.DataFrame({"Photon energy": list(counts_vs_photon_energy.keys()),
                                                "counts": list(counts_vs_photon_energy.values())})
                        print(df_terx["Photon energy"][0])
                        position_of_peak_terx = []

                        for i in range(len(df_terx["counts"]) // 100):
                            # problem jest taki że funkcja find_peak znajduje peaki po warunkach globalnych np. wyższe
                            # niż 20, ale chciałabym mieć średnie tło, tne na kawałki i w tych kawałkach określam liczbę
                            df_terx_cut = df_terx["counts"][i * 100:(i + 1) * 100]
                            average_terx = df_terx_cut.mean(axis=0)
                            for peak in find_peaks(df_terx_cut, height=(average_terx + 10 * average_terx ** .5))[0]:
                                position_of_peak_terx.append(float(df_terx["Photon energy"][i * 100 + peak]))

                        print(position_of_peak_terx, len(position_of_peak_terx))
                        energy_spx = []
                        counts_spx = []

                        for spx_line in my_spx.readlines():
                            if spx_line[0:4] == "Real":
                                measure_time = float(
                                    spx_line.strip("Real time:").strip()) / 3600000  # czas pomiaru w godzinach
                                print("measure_time", measure_time)
                            try:
                                float(spx_line.strip()[0:3])
                            except ValueError:
                                continue
                            spx_line_splited = [float(x) * 1000 for x in spx_line.strip().split()]

                            if E0_first * 1000 - 300 < spx_line_splited[0] < E0:
                                # chcę szukać peaków w pobliżu tych które znalazłam dla TERX'a

                                energy_spx.append(spx_line_splited[0])
                                counts_spx.append(spx_line_splited[1])

                        df_spx = pd.DataFrame({"Photon energy (SPX)": energy_spx, "Counts (SPX)": counts_spx})

                        position_of_peak_spx = []
                        for i in range(len(df_spx["Counts (SPX)"]) // 100):
                            average_spx = df_spx["Counts (SPX)"][i * 100:(i + 1) * 100].mean(axis=0)
                            for peak in find_peaks(df_spx["Counts (SPX)"][i * 100:(i + 1) * 100], height=(average_spx + 40 * (average_spx ** .5)))[0]:
                                position_of_peak_spx.append(df_spx["Photon energy (SPX)"][i*100+peak])

                        print(position_of_peak_spx, len(position_of_peak_spx))

                        #
                        #     try:
                        #         float(spx_line.split("   ").strip())
                        #     except ValueError:
                        #         continue
                        #     else:
                        # E_spx, count_spx = [float(x) for x in spx_line[0].strip("\n").split("\t")]
                        # print(E_spx)

                        # with open('testowy_rezultat.json', 'w') as t:
                        #     t.write(json.dumps(counts_vs_photon_energy))

    print('Calibration performed')



if __name__ == '__main__':
    calibration()

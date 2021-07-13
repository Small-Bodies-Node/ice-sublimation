import csv
from fastrot import *

def replicate_results():
    a_v_in = 0.05
    a_ir_in = 0.00
    inc_in = 90

    rh_list = [0.10, 0.13, 0.16, 0.20, 0.25, 0.32, 0.40, 0.50, 0.63,
               0.79, 1.00, 1.26, 1.58, 2.00, 2.51, 3.16, 3.98, 5.01,
               6.31, 7.94, 10.00, 12.59, 15.85, 19.95, 25.12]

    for species in speciesList:
        outfile = '_'.join(['fastrot', species.lower(), 'allr.csv'])
        with open('/'.join(['python-data', outfile]), "w") as file:
            wtr = csv.writer(file)
            wtr.writerow(('Species', 'Obl', 'r_H', 'log(r)', 'A_vis', 'A_ir', '<Z>', 'Log(<Z>)'))
            for rh_in in rh_list:
                nr = run_model(species, a_v_in, a_ir_in, rh_in, inc_in, )
                wtr.writerow((nr[0],
                              '{:.1f}'.format(nr[1]),
                              '{:.2f}'.format(nr[2]),
                              '{:.3f}'.format(nr[3]),
                              '{:.2f}'.format(nr[4]),
                              '{:.2f}'.format(nr[5]),
                              '{:.3e}'.format(nr[6]),
                              '{:.3f}'.format(nr[7]))
                             )

    return


if __name__ == '__main__':
    replicate_results()
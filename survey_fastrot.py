"""
    Description
    -----------
    Runs fastrot.py over a desired parameter space.

    The calculation results are stored in the files `results/output.json` and `results/output.csv`.
"""
import os
import csv


def survey_fastrot(species_list, Av_list, Air_list, rh_list, obl_list):
    import fastrot
    from itertools import product
    from json import dump

    search_space = []
    arguments = locals()
    for param in [species_list, Av_list, Air_list, rh_list, obl_list]:
        search_space.append(param if isinstance(param, list) else [param])
    results = []
    for inputs in product(*search_space):
        results.append(fastrot.run_model(*inputs))

    output_json = {'results': results}
    json_path = os.path.join(os.getcwd(), 'results', 'output.json')
    with open(json_path, 'w') as json_file:
        dump(output_json, json_file)

    fieldnames = results[0].keys()
    csv_path = os.path.join(os.getcwd(), 'results', 'output.csv')
    with open(csv_path, 'w') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    return output_json


description = "Use this program to iterate `fastrot.py` over a desired parameter space."

if __name__ == '__main__':
    speciesList = ['H2O', 'H2O-CH4', 'CO2', 'CO']

    import argparse

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--species_list', nargs='+', metavar='species', type=str, required=True,
                        help="Desired ice species to be considered. \n"
                             "The valid inputs are: " + ', '.join(speciesList)),
    parser.add_argument('--Av_list', nargs='+', metavar='visual albedo', type=float, required=True, )
    parser.add_argument('--Air_list', nargs='+', metavar='infrared albedo', type=float, required=True, )
    parser.add_argument('--rh_list', nargs='+', metavar='heliocentric_distance', type=float, required=True, )
    parser.add_argument('--obl_list', nargs='+', metavar='obliquity', type=float, required=True, )

    try:
        args = parser.parse_args()
        survey_fastrot(args.species_list, args.Av_list, args.Air_list, args.rh_list, args.obl_list)
    except Exception as e:
        print(e)

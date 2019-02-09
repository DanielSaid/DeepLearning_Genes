#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import sys
import pickle


class CompoundLabelParser:
    __COMPOUND_ABBREV = 'Abbreviation'
    __GENOTOXIC = 'GTX'
    __CARCINOGENIC = 'C'

    __GENOTOXIC_VITRO = 'GTX(vitro)'
    __GENOTOXIC_VIVO = 'GTX(vivo)'

    __GENOTOXIC_TR = 'GTX'
    __GENOTOXIC_F = 'NGTX'

    __CARCINOGENIC_TR = 'C'
    __CARCINOGENIC_F = 'NC'

    __UNCLASSIFIED = 'UC'

    def __init__(self):
        pass

    def parse(self, data_filename: str):
        labels_vitro = {}
        labels_vivo = {}
        with open(data_filename, errors='ignore') as data_stream:
            reader = csv.reader(data_stream, delimiter=';', quotechar='"')
            header = next(reader)
            idx_abbrev = header.index(CompoundLabelParser.__COMPOUND_ABBREV)
            idx_genotoxic = header.index(CompoundLabelParser.__GENOTOXIC)
            idx_carcinogenic = header.index(CompoundLabelParser.__CARCINOGENIC)
            ix_gtx_vt = header.index(CompoundLabelParser.__GENOTOXIC_VITRO)
            ix_gtx_vv = header.index(CompoundLabelParser.__GENOTOXIC_VIVO)

            for compound in reader:
                labels_vitro[compound[idx_abbrev]] = {}
                labels_vivo[compound[idx_abbrev]] = {}
                if compound[idx_genotoxic] != CompoundLabelParser.__UNCLASSIFIED:
                    labels_vitro[compound[idx_abbrev]][CompoundLabelParser.__GENOTOXIC] \
                        = compound[ix_gtx_vt] == CompoundLabelParser.__GENOTOXIC_TR
                    labels_vivo[compound[idx_abbrev]][CompoundLabelParser.__GENOTOXIC] \
                        = compound[ix_gtx_vv] == CompoundLabelParser.__GENOTOXIC_TR

                if compound[idx_carcinogenic] != CompoundLabelParser.__UNCLASSIFIED:
                    labels_vitro[compound[idx_abbrev]][CompoundLabelParser.__CARCINOGENIC] \
                        = compound[idx_carcinogenic] == CompoundLabelParser.__CARCINOGENIC_TR
                    labels_vivo[compound[idx_abbrev]][CompoundLabelParser.__CARCINOGENIC] \
                        = compound[idx_carcinogenic] == CompoundLabelParser.__CARCINOGENIC_TR
        return labels_vitro, labels_vivo


def main():
    if len(sys.argv) != 2:
        print("Please specify exactly 1 argument: [datafile]")
        sys.exit(1)

    parser = CompoundLabelParser()
    vitro, vivo = parser.parse(*sys.argv[1:])
    with open("labels_vitro.p", "wb") as outfile:
        pickle.dump(vitro, outfile)

    with open("labels_vivo.p", "wb") as outfile:
        pickle.dump(vivo, outfile)


if __name__ == '__main__':
    main()

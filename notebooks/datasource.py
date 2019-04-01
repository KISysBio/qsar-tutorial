"""QSAR modelling Data Sources

This module supports representation of functional bioactivies

- Directly from ChEMBL API

- From a list of Canonical SMILES, along with the measured activity

"""

import numpy as np
import pandas as pd

from chembl_webresource_client.new_client import new_client


class Preprocessing:
    """Preprocess data obtained via a DataSource (.csv, .xls file or ChEMBL)

    The bioactivity data is filtered, duplicated and invalid entries are handled and
      molecular descriptors are calculated.
    """

    def __init__(self, activity_column='median_activity_value'):
        self.activity_column = activity_column

    def _mark_to_remove(self, data, activity_column):
        """Mark data for removal if standard deviation of duplicated entries is above 1"""

        activity_value = data[activity_column].astype(float)
        std_deviation = np.std(activity_value)
        number_compounds = data.shape[0]

        # It is duplicated
        result_df = pd.Series(
            {'number_compounds': number_compounds,
             'standard_deviation': std_deviation,
             self.activity_column: np.median(activity_value),
             'duplicated': number_compounds > 1,
             'mark_to_remove': std_deviation > 1})
        return result_df

    def _default_filter_chembl(self, bioactivities_df):
        """Default filter for data downloaded from ChEMBL"""

        # Relation measured must be of type equality
        valid_rows = bioactivities_df['relation'] == '='
        bioactivities_df = bioactivities_df[valid_rows]

        # Remove data that is marked as potential invalid
        valid_rows = bioactivities_df['data_validity_comment'].isnull()
        bioactivities_df = bioactivities_df[valid_rows]

        # There must be a valid pCHEMBL value
        valid_rows = ~bioactivities_df['pchembl_value'].isnull()
        bioactivities_df = bioactivities_df[valid_rows]
        return bioactivities_df

    def clean_dataset(self, data_source):

        bioactivities_df = data_source.bioactivities_df.copy()

        if type(data_source) is ChEMBLApiDataSource:
            # Perform initial filter on bioactivities
            bioactivities_df = self._default_filter_chembl(bioactivities_df)
        else:
            raise ValueError("Data source not recognized: %s" % type(data_source))

        # Mark duplicated entries for removeal
        grouped_dataset = bioactivities_df.groupby(data_source.compound_id_column)
        grouped_dataset = grouped_dataset.apply(lambda x: self._mark_to_remove(x, data_source.activity_column))

        # The resulting dataset, clean_df, does not contain duplicated entries
        merged_df = pd.merge(bioactivities_df, grouped_dataset.reset_index())
        merged_df = merged_df[~merged_df['mark_to_remove']]
        clean_df = merged_df.groupby(data_source.compound_id_column).head(1)

        return clean_df


class DataSource():
    """
    Represent a data source

    Chemical data might be obtained from different sources:

        - ChEMBL API (class ChEMBLApiDataSource)
        - A CSV file downloaded from ChEMBL (class ChEMBLCsvDataSource)
        - A CSV file produced manually or exported from elsewhere containing two columns:
            'canonical_smiles' and 'activity' (class CsvDataSource)
          The first represent the Canonical Smiles of a compound while the second correspond
            to the acitivity of the compounds.

        DataSource parses the data and creates an object of the class modSAR.qsar.QSARDataset
          using the function `get_qsar_dataset`.

    """

    def get_qsar_dataset(self):
        raise NotImplementedError()


class ChEMBLApiDataSource(DataSource):
    """
        Retrieve active compounds for a target using ChEMBL API

        Example of use:

            chembl_data_source = ChEMBLApiDataSource(target_id='CHEMBL202', standard_types=['IC50'])
            chembl202_dataset = chembl_data_source.get_qsar_dataset()
    """

    def __init__(self, target_id, standard_types):
        """
        Stores a DataFrame containing the bioactivities listed on ChEMBL for specified target.
        Bioactivities can later be converted to a QSARDataset using function `get_qsar_dataset`

        Args:
            target_id (str): ChEMBL target id
            standard_types (str or list): e.g.: IC50, Ki, etc.

        """

        self.target_id = target_id
        if type(standard_types) is str:
            self.standard_types = [standard_types]
        else:
            self.standard_types = standard_types

        activity = new_client.activity
        result = activity.filter(target_chembl_id=target_id, assay_type__iregex='(B|F)')

        def get_compound_df(idx, compound_dict):
            """Filter compounds by the standard_types informed, returning a DataFrame"""

            if self.standard_types is None or compound_dict['standard_type'] is None:
                is_valid_std_type = False
            else:
                compound_std_type = compound_dict['standard_type'].lower()
                is_valid_std_type = any([std_type.lower() in compound_std_type
                                         for std_type in self.standard_types])

            if is_valid_std_type:
                # Drop unused columns
                compound_dict.pop('activity_properties', None)
                # Capture Ligand Efficiency
                lig_efficiency = compound_dict.pop('ligand_efficiency', None)

                compound_df = pd.DataFrame(compound_dict, index=[0])
                if lig_efficiency is not None:
                    lig_efficiency_df = pd.DataFrame(lig_efficiency, index=[0])
                    lig_efficiency_df.columns = ['ligand_efficiency_%s' % col
                                                 for col in lig_efficiency_df.columns]
                    compound_df = pd.concat([compound_df, lig_efficiency_df], axis=1)
            else:
                compound_df = pd.DataFrame(columns=compound_dict.keys())

            print_progress_bar(idx + 1, self.number_retrieved_compounds,
                               prefix='Progress:', suffix='Complete', length=50)
            return compound_df

        self.number_retrieved_compounds = len(result)
        print_progress_bar(0, self.number_retrieved_compounds,
                           prefix='Progress:', suffix='Complete', length=50)
        bioactivities_df = [get_compound_df(idx, compound_dict)
                            for idx, compound_dict in enumerate(result)]
        bioactivities_df = pd.concat(bioactivities_df, sort=False).reset_index(drop=True)

        self.bioactivities_df = bioactivities_df

        self.smiles_column = 'canonical_smiles'
        self.compound_id_column = 'parent_molecule_chembl_id'
        self.activity_column = 'pchembl_value'

    def save_bioactivities(self, xls_filename):
        self.bioactivities_df.to_excel(xls_filename, index=False)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        repr_str = \
            "ChEMBLApiDataSource object" \
            + "\n  target_id: " + self.target_id \
            + "\n  bioactivities: %d" % self.bioactivities_df.shape[0] \
            + "\n  standard_types: " + str(self.bioactivities_df['standard_type'].unique())
        return repr_str


# Print iterations progress
def print_progress_bar(iteration, total, prefix='', suffix='',
                       decimals=1, length=100, fill='█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end='\r')
    # Print New Line on Complete
    if iteration == total:
        print()

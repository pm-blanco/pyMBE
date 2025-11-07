import pandas as pd
import json
import re
import numpy as np
import logging
import warnings

class _DFManagement:

    class _NumpyEncoder(json.JSONEncoder):
        """
        Custom JSON encoder that converts NumPy arrays to Python lists
        and NumPy scalars to Python scalars.
        """
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.generic):
                return obj.item()
            return super().default(obj)

    @classmethod
    def _add_bond_in_df(cls, df, particle_id1, particle_id2, use_default_bond=False):
        """
        Adds a bond entry on the `pymbe.df` storing the particle_ids of the two bonded particles.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            particle_id1(`int`): particle_id of the type of the first particle type of the bonded particles
            particle_id2(`int`): particle_id of the type of the second particle type of the bonded particles
            use_default_bond(`bool`, optional): Controls if a bond of type `default` is used to bond particle whose bond types are not defined in `pmb.df`. Defaults to False.

        Returns:
            df(`DataFrame`): dataframe with pyMBE information with the new bond added.
            index(`int`): Row index where the bond information has been added in pmb.df.
        """
        particle_name1 = df.loc[(df['particle_id']==particle_id1) & (df['pmb_type']=="particle")].name.values[0]
        particle_name2 = df.loc[(df['particle_id']==particle_id2) & (df['pmb_type']=="particle")].name.values[0]

        bond_key = cls._find_bond_key(df = df, 
                                      particle_name1 = particle_name1,
                                      particle_name2 = particle_name2,
                                      use_default_bond = use_default_bond)
        if not bond_key:
            return None
        df = cls._copy_df_entry(df = df,
                                name = bond_key,
                                column_name = 'particle_id2',
                                number_of_copies = 1)
        indexs = np.where(df['name'] == bond_key)
        index_list = list(indexs[0])
        used_bond_df = df.loc[df['particle_id2'].notnull()]
        #without this drop the program crashes when dropping duplicates because the 'bond' column is a dict
        used_bond_df = used_bond_df.drop([('bond_object','')],axis =1 )
        used_bond_index = used_bond_df.index.to_list()
        if not index_list:
            return None
        for index in index_list:
            if index not in used_bond_index:
                cls._clean_df_row(df = df,
                                  index = int(index))
                df.at[index,'particle_id'] = particle_id1
                df.at[index,'particle_id2'] = particle_id2
                break
        return df, index

    @classmethod
    def _add_value_to_df(cls, df, index,key,new_value, non_standard_value=False, overwrite=False):
        """
        Adds a value to a cell in the `pmb.df` DataFrame.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            index(`int`): index of the row to add the value to.
            key(`str`): the column label to add the value to.
            non_standard_value(`bool`, optional): Switch to enable insertion of non-standard values, such as `dict` objects. Defaults to False.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False.
        """

        token = "#protected:"

        def protect(obj):
            if non_standard_value:
                return token + json.dumps(obj, cls=cls._NumpyEncoder)
            return obj

        def deprotect(obj):
            if non_standard_value and isinstance(obj, str) and obj.startswith(token):
                return json.loads(obj.removeprefix(token))
            return obj

        # Make sure index is a scalar integer value
        index = int(index)
        assert isinstance(index, int), '`index` should be a scalar integer value.'
        idx = pd.IndexSlice
        if cls._check_if_df_cell_has_a_value(df=df, index=index, key=key):
            old_value = df.loc[index,idx[key]]
            if not pd.Series([protect(old_value)]).equals(pd.Series([protect(new_value)])):
                name= df.loc[index,('name','')]
                pmb_type= df.loc[index,('pmb_type','')]
                logging.debug(f"You are attempting to redefine the properties of {name} of pmb_type {pmb_type}")    
                if overwrite:
                    logging.info(f'Overwritting the value of the entry `{key}`: old_value = {old_value} new_value = {new_value}')
                if not overwrite:
                    logging.debug(f"pyMBE has preserved of the entry `{key}`: old_value = {old_value}. If you want to overwrite it with new_value = {new_value}, activate the switch overwrite = True ")
                    return

        df.loc[index,idx[key]] = protect(new_value)
        if non_standard_value:
            df[key] = df[key].apply(deprotect)
        return
    
    @classmethod
    def _assign_molecule_id(cls, df, molecule_index):
        """
        Assigns the `molecule_id` of the pmb object given by `pmb_type`
        
        Args:
            molecule_index(`int`): index of the current `pmb_object_type` to assign the `molecule_id`
        Returns:
            molecule_id(`int`): Id of the molecule
        """
        cls._clean_df_row(df = df,
                          index = int(molecule_index))

        if df['molecule_id'].isnull().values.all():
            molecule_id = 0        
        else:
            molecule_id = df['molecule_id'].max() +1
        cls._add_value_to_df(df = df,
                             key = ('molecule_id',''),
                             index = int(molecule_index),
                             new_value = molecule_id)
        return molecule_id

    @staticmethod
    def _check_if_df_cell_has_a_value(df, index, key):
        """
        Checks if a cell in the `pmb.df` at the specified index and column has a value.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            index(`int`): Index of the row to check.
            key(`str`): Column label to check.

        Returns:
            `bool`: `True` if the cell has a value, `False` otherwise.
        """
        idx = pd.IndexSlice
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return not pd.isna(df.loc[index, idx[key]])

    @staticmethod
    def _check_if_name_is_defined_in_df(name, df):
        """
        Checks if `name` is defined in `pmb.df`.

        Args:
            name(`str`): label to check if defined in `pmb.df`.
            df(`DataFrame`): dataframe with pyMBE information.

        Returns:
            `bool`: `True` for success, `False` otherwise.
        """
        return name in df['name'].unique()

    @staticmethod
    def _check_if_multiple_pmb_types_for_name(name, pmb_type_to_be_defined, df):
        """
        Checks if `name` is defined in `pmb.df` with multiple pmb_types.

        Args:
            name(`str`): label to check if defined in `pmb.df`.
            pmb_type_to_be_defined(`str`): pmb object type corresponding to `name`.
            df(`DataFrame`): dataframe with pyMBE information.

        Returns:
            `bool`: `True` for success, `False` otherwise.
        """
        if name in df['name'].unique():
            current_object_type = df[df['name']==name].pmb_type.values[0]
            if current_object_type != pmb_type_to_be_defined:
                raise ValueError (f"The name {name} is already defined in the df with a pmb_type = {current_object_type}, pymMBE does not support objects with the same name but different pmb_types")

    @classmethod
    def _clean_df_row(cls, df, index, columns_keys_to_clean=("particle_id", "particle_id2", "residue_id", "molecule_id")):
        """
        Cleans the columns of `pmb.df` in `columns_keys_to_clean` of the row with index `index` by assigning them a pd.NA value.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            index(`int`): Index of the row to clean.
            columns_keys_to_clean(`list` of `str`, optional): List with the column keys to be cleaned. Defaults to [`particle_id`, `particle_id2`, `residue_id`, `molecule_id`].
        """   
        for column_key in columns_keys_to_clean:
            cls._add_value_to_df(df = df,
                                 key = (column_key,''),
                                 index = index,
                                 new_value = pd.NA)
        df.fillna(pd.NA, 
                  inplace = True)

    @staticmethod
    def _clean_ids_in_df_row(df, row):
        """
        Cleans particle, residue and molecules ids in `row`.
        If there are other repeated entries for the same name, drops the row.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            row(pd.DataFrame): A row from the DataFrame to clean.

        Returns:
            df(`DataFrame`): dataframe with pyMBE information with cleaned ids in `row    
        """
        columns_to_clean = ['particle_id',
                            'particle_id2', 
                            'residue_id', 
                            'molecule_id']
        if len(df.loc[df['name'] == row['name'].values[0]]) > 1:
            df = df.drop(row.index).reset_index(drop=True)
            
        else:
            for column_name in columns_to_clean:
                df.loc[row.index, column_name] = pd.NA
        return df

    @staticmethod
    def _copy_df_entry(df, name, column_name, number_of_copies):
        '''
        Creates 'number_of_copies' of a given 'name' in `pymbe.df`.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            name(`str`): Label of the particle/residue/molecule type to be created. `name` must be defined in `pmb.df`
            column_name(`str`): Column name to use as a filter. 
            number_of_copies(`int`): number of copies of `name` to be created.

        Returns:
            df(`DataFrame`): dataframe with pyMBE information with the new copies of `name` added.

        Note:
            - Currently, column_name only supports "particle_id", "particle_id2", "residue_id" and "molecule_id" 
        '''
        valid_column_names=["particle_id", "residue_id", "molecule_id", "particle_id2" ]
        if column_name not in valid_column_names:
            raise ValueError(f"{column_name} is not a valid column_name, currently only the following are supported: {valid_column_names}")
        df_by_name = df.loc[df.name == name]
        if number_of_copies != 1:                           
            df_by_name_repeated = pd.concat ([df_by_name]*(number_of_copies-1), ignore_index=True)
            # Concatenate the new particle rows to  `df`
            df = pd.concat ([df,df_by_name_repeated], ignore_index=True)
        else:
            if not df_by_name[column_name].isnull().values.any():     
                df_by_name = df_by_name[df_by_name.index == df_by_name.index.min()] 
                df_by_name_repeated = pd.concat ([df_by_name]*(number_of_copies), ignore_index=True)
                df_by_name_repeated[column_name] = pd.NA
                df = pd.concat ([df,df_by_name_repeated], ignore_index=True)
        return df

    @staticmethod
    def _create_variable_with_units(variable, units_registry):
        """
        Returns a pint object with the value and units defined in `variable`.

        Args:
            variable(`dict` or `str`): {'value': value, 'units': units}
            units_registry(`pint.UnitRegistry`): pyMBE UnitRegistry object.

        Returns:
            variable_with_units(`obj`): variable with units using the pyMBE UnitRegistry.
        """        
        if isinstance(variable, dict):
            value=variable.pop('value')
            units=variable.pop('units')
        elif isinstance(variable, str):
            value = float(re.split(r'\s+', variable)[0])
            units = re.split(r'\s+', variable)[1]
        variable_with_units = value * units_registry(units)
        return variable_with_units

    @classmethod
    def _convert_columns_to_original_format(cls,df,units_registry):
        """
        Converts the columns of the Dataframe to the original format in pyMBE.
        
        Args:
            df(`DataFrame`): dataframe with pyMBE information as a string
            units_registry(`pint.UnitRegistry`): pyMBE UnitRegistry object.  
        
        """

        columns_dtype_int = ['particle_id','particle_id2', 'residue_id','molecule_id', ('state_one','es_type'),('state_two','es_type'),('state_one','z'),('state_two','z') ]  

        columns_with_units = ['sigma', 'epsilon', 'cutoff', 'offset']

        columns_with_list_or_dict = ['residue_list','side_chains', 'parameters_of_the_potential','sequence', 'chain_map', 'node_map']

        for column_name in columns_dtype_int:
            df[column_name] = df[column_name].astype(pd.Int64Dtype())
            
        for column_name in columns_with_list_or_dict:
            if df[column_name].isnull().all():
                df[column_name] = df[column_name].astype(object)
            else:
                df[column_name] = df[column_name].apply(lambda x: json.loads(x) if pd.notnull(x) else x)

        for column_name in columns_with_units:
            df[column_name] = df[column_name].apply(lambda x: cls._create_variable_with_units(x, units_registry) if pd.notnull(x) else x)

        df['bond_object'] = df['bond_object'].apply(lambda x: cls._convert_str_to_bond_object(x) if pd.notnull(x) else x)
        df["l0"] = df["l0"].astype(object)
        df["pka"] = df["pka"].astype(object)

    @staticmethod
    def _convert_str_to_bond_object(bond_str):
        """
        Convert a row read as a `str` to the corresponding ESPResSo bond object. 

        Args:
            bond_str(`str`): string with the information of a bond object.

        Returns:
            bond_object(`obj`): ESPResSo bond object.

        Note:
            Current supported bonds are: HarmonicBond and FeneBond
        """
        import espressomd.interactions

        supported_bonds = ['HarmonicBond', 'FeneBond']
        m = re.search(r'^([A-Za-z0-9_]+)\((\{.+\})\)$', bond_str)
        if m is None:
            raise ValueError(f'Cannot parse bond "{bond_str}"')
        bond = m.group(1)
        if bond not in supported_bonds:
            raise NotImplementedError(f"Bond type '{bond}' currently not implemented in pyMBE, accepted types are {supported_bonds}")
        params = json.loads(m.group(2))
        bond_id = params.pop("bond_id")
        bond_object = getattr(espressomd.interactions, bond)(**params)
        bond_object._bond_id = bond_id
        return bond_object

    @staticmethod
    def delete_entries_in_df(df, entry_name):
        """
        Deletes entries with name `entry_name` from the DataFrame if it exists.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            entry_name (`str`): The name of the entry in the dataframe to delete.

        Returns:
            df(`DataFrame`): dataframe with pyMBE information with the entry deleted.
        """
        if entry_name in df["name"].values:
            df = df[df["name"] != entry_name].reset_index(drop=True)
        return df

    @staticmethod
    def _find_bond_key(df, particle_name1, particle_name2, use_default_bond=False):
        """
        Searches for the `name` of the bond between `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        
        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            particle_name1(`str`): label of the type of the first particle type of the bonded particles.
            particle_name2(`str`): label of the type of the second particle type of the bonded particles.
            use_default_bond(`bool`, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to 'False'. 

        Returns:
            bond_key (str): `name` of the bond between `particle_name1` and `particle_name2` if a matching bond exists

        Note:
            - If `use_default_bond`=`True`, it returns "default" if no key is found.
        """
        bond_keys = [f'{particle_name1}-{particle_name2}', f'{particle_name2}-{particle_name1}']
        bond_defined=False
        for bond_key in bond_keys:
            if bond_key in df["name"].values:
                bond_defined=True
                correct_key=bond_key
                break
        if bond_defined:
            return correct_key
        elif use_default_bond:
            return 'default'
        else:
            return None

    @staticmethod
    def _setup_df():
        """
        Sets up the pyMBE's dataframe `pymbe.df`.

        Returns:
            columns_names(`obj`): pandas multiindex object with the column names of the pyMBE's dataframe
        """
        
        columns_dtypes = {
            'name': {
                '': str},
            'pmb_type': {
                '': str},
            'particle_id': {
                '': pd.Int64Dtype()},
            'particle_id2':  {
                '': pd.Int64Dtype()},
            'residue_id':  {
                '': pd.Int64Dtype()},
            'molecule_id':  {
                '': pd.Int64Dtype()},
            'acidity':  {
                '': str},
            'pka':  {
                '': object},
            'central_bead':  {
                '': object},
            'side_chains': {
                '': object},
            'residue_list': {
                '': object},
            'model': {
                '': str},
            'sigma': {
                '': object},
            'cutoff': {
                '': object},
            'offset': {
                '': object},
            'epsilon': {
                '': object},
            'state_one': {
                'label': str,
                'es_type': pd.Int64Dtype(),
                'z': pd.Int64Dtype()},
            'state_two': {
                'label': str,
                'es_type': pd.Int64Dtype(),
                'z': pd.Int64Dtype()},
            'sequence': {
                '': object},
            'bond_object': {
                '': object},
            'parameters_of_the_potential':{
                '': object},
            'l0': {
                '': float},
            'node_map':{
                '':object},
            'chain_map':{
                '':object}}
        
        df = pd.DataFrame(columns=pd.MultiIndex.from_tuples([(col_main, col_sub) for col_main, sub_cols in columns_dtypes.items() for col_sub in sub_cols.keys()]))
        
        for level1, sub_dtypes in columns_dtypes.items():
            for level2, dtype in sub_dtypes.items():
                df[level1, level2] = df[level1, level2].astype(dtype)

        columns_names = pd.MultiIndex.from_frame(df)
        columns_names = columns_names.names
                
        return df
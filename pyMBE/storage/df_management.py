import pandas as pd
import json
import re

class _DFManagement:

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
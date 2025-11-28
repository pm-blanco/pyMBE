#
# Copyright (C) 2024 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import pandas as pd
import numpy as np

def add_data_to_df(df, data_dict, index):
    """
    Adds the data in data_dict in the df pandas dataframe
    
    Args:
        df(`Pandas.Dataframe`): pandas dataframe
        data_dict(`dict`): dictionary with the data to be added
        index(`lst`): index of the df where the data should be added
    Returns:
        updated_df(`Pandas.Dataframe`): pandas dataframe updated with data_dict

    """
    if df.empty:
        # if the dataframe is empty, load the first set of data
        updated_df = pd.DataFrame(data_dict,
                          index=index)
    else:
        # if the dataframe has data, concatenate the new data
        updated_df = pd.concat([df, pd.DataFrame(data_dict,
                         index=index)])
    return updated_df

def analyze_time_series(path_to_datafolder, filename_extension= ".csv", minus_separator = False, ignore_files=None):
    """
    Analyzes all time series stored in `path_to_datafolder` using the block binning method.

    Args:
        path_to_datafolder(`str`): path to the folder with the files with the time series
        filename_extension(`str`): extension of the file. Defaults to ".csv"
        minus_separator(`bool`): switch to enable the minus as a separator. Defaults to False.
        ignore_files(`lst`): list of filenames to be ignored for the binning analysis.

    Returns:
        data(`Pandas.Dataframe`): Dataframe with the time averages of all the time series in the datafolder.

    Note:
        - For more information about `minus_separator`, please check the documentation of `get_params_from_file_name`

    """
    data=pd.DataFrame()
    if ignore_files is None:
        ignore_files=[]
    with os.scandir(path_to_datafolder) as subdirectory:
        # Gather all data
        for subitem in subdirectory:
            if subitem.is_file():
                ignore_file=False
                for file in ignore_files:
                    if set(file.split()) == set(subitem.name.split()):
                        ignore_file=True
                if ignore_file:
                    continue
                if filename_extension in subitem.name:
                    # Get parameters from the file name
                    data_dict=get_params_from_file_name(file_name=subitem.name,
                                                        filename_extension=filename_extension,
                                                        minus_separator=minus_separator)
                    # Get the observables for binning analysis
                    time_series_data=pd.read_csv(f"{path_to_datafolder}/{subitem.name}")
                    analyzed_data=block_analyze(full_data=time_series_data)
                    value_list=[]
                    index_list=[]
                    for key in data_dict.keys():
                        value_list.append(data_dict[key])
                        index_list.append((key,"value"))
                    analyzed_data = pd.concat([pd.Series(value_list, index=index_list), analyzed_data])
                    data = add_data_to_df(df=data,
                                        data_dict=analyzed_data.to_dict(),
                                        index=[len(data)])   
    return data

def block_analyze(full_data, n_blocks=16, time_col = "time", equil=0.1,  columns_to_analyze = "all", verbose = False, dt=None):
    """
    Analyzes the data in `full_data` using a binning procedure.

    Args:
        full_data(`Pandas.Dataframe`): dataframe with the observables time series
        n_blocks(`int`): number of blocks used in the binning procedure.
        time_col(`str`): column name where the time is stored in `full_data`. Defaults to `time`.
        equil(`float`,opt): fraction of the data discarded as part of the equilibration. Defaults to 0.1.
        columns_to_analyze(`list`): array of column names to be analyzed. Defaults to "all".
        verbose(`bool`): switch to activate/deactivate printing the block size. Defaults to False.
        dt(`float`): time step in which the data was stored. If no value is provided, is calculated from the time series. Defaults to None.

    Returns:
        `result`(`Pandas.Series`): Series the mean (mean), statistical error (err_mean), number of effective samples (n_eff) and correlation time (tau_int) of each observable.
    """
    if dt is None:
        dt, n_warnings = get_dt(full_data, time_col, verbose = verbose) # check that the data was stored with the same time interval dt
        if n_warnings > 0 and verbose:
            print(f"Warning: looks like a repeated time value was encountered {n_warnings} times")

    drop_rows = int(full_data.shape[0]*equil) # calculate how many rows should be dropped as equilibration
    # drop the rows that will be discarded as equilibration
    data = full_data.drop(range(0,drop_rows))
    # drop the columns step, time and MC sweep
    if time_col in data.columns :
        data = data.drop(columns = time_col)
    else:
        raise ValueError(f"could not find the time column {time_col} in the data")
    if columns_to_analyze != "all":
        for column_name in data.columns:
            if column_name not in columns_to_analyze:
                data = data.drop(columns=column_name)
    # first, do the required operations on the remaining data
    n_samples = data.shape[0] # number of samples to be analyzed
    block_size = n_samples/n_blocks # mean block size
    mean = data.mean() # calculate the mean values of each column
    var_all = data.var() # calculate the variance of each column
    if verbose:
        print(f"n_blocks b = {n_blocks},  block_size k = {block_size}")

    # calculate the mean per each block
    blocks = split_dataframe(df=data,
                             n_blocks=n_blocks)

    block_means = [] # empty list that we will use to store the means per each block
    for block in blocks:
        block_mean = block.mean() # mean values of each column within a given block
        block_means.append(block_mean)
    block_means = pd.concat(block_means, axis=1).transpose()

    # perform calculations using averages or individual data blocks
    var_blocks = (n_blocks)/(n_blocks-1)*block_means.var() # variance of the block averages = variance of the mean
    err_mean = np.sqrt(var_blocks/n_blocks) # standard error of the mean by eq.(37) of Janke
    tau_int = dt*block_size * var_blocks / var_all /2.# eq.(38) of Janke
    n_eff = n_samples / (2*tau_int/dt) # effective number of samples in the whole simulation using eq.(28) of Janke

    # concatenate the observables and atribute a key for each (unique index)
    result = pd.concat( [ mean, err_mean, n_eff, tau_int], keys= [ "mean", "err_mean", "n_eff", "tau_int" ], join="inner")
    
    return result

import re

def built_output_name(input_dict):
    """
    Builts the output name for a given set of input parameters.
    
    Args:
        input_dict(`dict`): dictionary with all terminal inputs.

    Returns:
        output_name(`str`): name used for the output files
    
    Note:
        - The standard formatting rule is parametername1_parametervalue1_parametername2_parametervalue2
    """
    output_name=""
    for label in input_dict:
        if isinstance(input_dict[label], (str, bool)):
            formatted_variable = f"{input_dict[label]}"
        else:  
            formatted_variable = f"{input_dict[label]:.3g}"
        
        formatted_variable = formatted_variable.replace(" ", "_")
        formatted_variable = re.sub(r"[^\w\-.]", "", formatted_variable)
        output_name += f"{label}_{formatted_variable}_"

    return output_name.rstrip("_")  # Remove trailing underscore

def get_dt(data, time_col = "time", relative_tolerance = 0.01, verbose = False):
    """ 
    Calculate the time difference between the lines of the time series

    Args:
        data(`obj`): PandasDataFrame containing the time series
        time_col(`str`): string denoting the name of the column which contains the time data
        relative_toleranace(`float`): threshold up to which an inconsistency between consecutive dt values is considered acceptable (see Notes below for more details)
        verbose(`bool`): Switch to activate/deactivate printing of warning messages.

    Returns:
        (`float`, `int`) 1. value of time difference between consecutive lines, 2. number of repeated time values

    Note:
        - `relative_tolerance` is used in two ways to check for irregularities in time intervals of between consecutive lines:
            - Sometimes, small inconsistencies in time intervals may occur due to rounding errors or because we run a short tuning of cells or Ewald parameters. These inconsistencies are considered negligible, as long as they do not exceed dt * relative_tolerance
            - Sometimes, consecutive values may have identical times, e.g. if only reactions are performed but there no MD evolution in between. In this case, dt is ill-defined. As of May 2024, we only produce warnings but do not raise an error in order to retain backward compatibility. We deprecate such use and plan to remove this compatibility in the future. If necessary, the user can always override get_dt() by providing a custom value on the input.
    """
    if time_col in data.columns.to_list():
        time = data[ time_col ]
    else:
        raise ValueError(f"Column \'{time_col}\' not found in columns: "+str( data.columns.to_list() ) )
    imax = data.shape[0]
    dt_init = time[1] - time[0]
    if dt_init < 1e-8:
        raise ValueError(f"The two first rows contain data samples at the same simulation time: time[0] = {time[0]} time[1] = {time[1]}. Post-processing of data with repeated time values is not supported because it breaks the estimation of the autocorrelation time.")
    warn_lines = []
    for i in range(1,imax):
        dt = time[i] - time[i-1]
        if np.abs( (dt_init - dt)/dt_init ) > relative_tolerance : # dt is different from before
            if np.abs( dt/dt_init ) < relative_tolerance : # looks like a repeated value
                warn_lines.append(f"Row {i} dt = {dt} = {time[i]} - {time[i-1]} not equal to dt_init = {dt_init}")
                dt = max(dt,dt_init)
            else:
                raise ValueError(f"\nget_dt: Row {i} dt = {dt} = {time[i]} - {time[i-1]} not equal to dt_init = {dt_init}\n")
    if verbose:
        print(warn_lines)
    return dt, len(warn_lines)

def get_params_from_file_name(file_name, minus_separator = False, filename_extension=".csv"):
    """
    Gets the parameters from `file_name`. By default, it assumes the standard formating rule  `file_name=obsname1_value1_obsname2_value2...`
    
    Args:
        file_name(`str`): name of the file
        minus_separator(`bool`): switch to enable the minus as a separator, see Notes. Defaults to False.
        filename_extension(`str`): extension of the file. Defaults to ".csv"
    
    Returns:
        params(`dict`): dictionary with the labels and values of the parameters.

    Note:
        - If `minus_separator = True`, then the data is parsed assuming the following formating rule  `file_name=obsname1-value1_obsname2-value2...`. 
        - This is functionality is kept for backwards compatibility.
    """
    file_name = os.path.basename(file_name)
    params = {}
    if filename_extension in file_name:
        system_name = file_name.replace(f"{filename_extension}", '')
    else:
        system_name = file_name
    if minus_separator:
        entries = system_name.split('_')
        for entry in entries:
            splitted_entry = entry.split('-')
            if splitted_entry[0] == 'N':
                params[splitted_entry[0]] = int(splitted_entry[-1])
            else:
                params[splitted_entry[0]] = splitted_entry[-1]           
    else:
        entries = system_name.split('_')
        if len(entries) % 2:
            raise ValueError("Wrong file name format. Need even number of entries separated by underscores, got: " + str(entries) + str( len(entries)) )
        for i in range( len(entries)//2 ):
            params[ entries[2*i] ] = entries[2*i+1]
    return params

def split_dataframe(df,n_blocks):
    """
    Splits a Pandas Dataframe in n_blocks of approximately the same size.

    Args:
        df(`Pandas.Dataframe`): Dataframe to be split
        n_blocks(`int`): Number of blocks

    Returns:
        (`lst` of `Pandas.Dataframe`): array with the splitted dataframes.

    Note:
        - For a `df` of length `l` that should be split into n_blocks, it returns l % n_blocks sub-arrays of size l//n_blocks + 1 and the rest of size l//n_blocks.
        - The behaviour of this function is the same as numpy.array_split for backwards compatibility, see [docs](https://numpy.org/doc/stable/reference/generated/numpy.array_split.html)

    """

    # Blocks of size 1 (s1) =  df.shape[0]//n_blocks+1

    n_blocks_s1=df.shape[0] % n_blocks
    block_size_s1=df.shape[0]//n_blocks+1
    blocks=split_dataframe_in_equal_blocks(df=df,
                                           start_row=0,
                                           end_row=n_blocks_s1*block_size_s1,
                                           block_size=block_size_s1)


    # Blocks of size 2 (s2) =  df.shape[0]//n_blocks
    block_size_s2=df.shape[0]//n_blocks
    blocks+=split_dataframe_in_equal_blocks(df=df,
                                           start_row=n_blocks_s1*block_size_s1,
                                           end_row=df.shape[0],
                                           block_size=block_size_s2)
    return blocks

def split_dataframe_in_equal_blocks(df, start_row, end_row, block_size):
    """
    Splits a Pandas dataframe in equally spaced blocks.

    Args:
        df(`Pandas.Dataframe`): Dataframe to be splitted.
        start_row(`int`): index of the first row.
        end_row(`int`): index of the last row.
        block_size(`int`): number of rows per block.

    Returns:
        (`lst`): array of dataframes of equal size.
    """
    return [df[row:row+block_size] for row in range(start_row,end_row,block_size)]









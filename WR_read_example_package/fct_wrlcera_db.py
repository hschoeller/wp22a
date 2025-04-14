#!/usr/bin/env python

# original code by Dominik Büeler dominik.bueeler@meteoswiss.ch 
# 20250102: Adaptations by christian.grams@gmx.de
#
###################
# Define function #
###################

# Function
def wrlcera(wr=['AT','ZO','ScTr','AR','EuBL','ScBL','GL','no'], start='19790101_00', end='20241231_21', tformat='string', setup='z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019', dataset='era5', basepath = './mv0561/ec.era5/eof/update_1950/latwgt/eof_AUTO_Z500_N161_-80E40E30N90N_year/norminput'):

    '''
    Description: Load individual life cycles of a specific weather regime for a
                 specific period from era-interim (for details about the weather regimes
                 contact Christian Grams, christian.grams@gmx.de)

    Input: wr     = list of weather regimes for which the life cycles should be loaded;
                    default = ['AT','ZO','ScTr','AR','EuBL','ScBL','GL','no']
           start  = start date as string in the format YYYYMMDD_HH that includes the first onset;
                    default = '19790101_00'
           end    = end date as string in the format YYYYMMDD_HH that includes the last decay;
                    default = '20241231_21'
           tformat = 'string' if dates should be returned in string format,
                     'dtime' if dates should be returned in datetime format;
                     default = 'string'
           setup = string defining the weather regime setup from which the data should
                   be read -> it is important to select this carefully!;
                   'z500anom_1979_2015_on_wrdef_10d_1.0_1979_2015' = z500 anomalies with respect
                   to climatology of 1979-2015 are projected on weather regime cluster means
                   based on 10d-low-pass-filtered z500 anomalies between 1979 and 2015
                   with an iwr threshold of 1.0 (erainterim),
                   'z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019' = z500 anomalies with respect
                   to climatology of 1979-2019 are projected on weather regime cluster means
                   based on 10d-low-pass-filtered z500 anomalies between 1979 and 2019
                   with an iwr threshold of 1.0 (era5),
           dataset = 'erainterim' or 'era5'; default = 'era5'
           basepath = directory to the raw data

    Output: Dictionary with the structure dict[key], with the keys being the different weather regimes;
            for each weather regime (except no regime), the following characteristic values
            are given for every individual life cycle:
                [0]:  'number':     life cycle id (only unique within the corresponding weather regime),
                [1]:  'onset':      date of onset,
                [2]:  'sat_start':  date of saturation start,
                [3]:  'mx':         date of maximum projection,
                [4]:  'sat_end':    date of saturation end,
                [5]:  'decay':      date of decay,
                [6]:  'dcfr':       of secondary importance; counterpart of onto
                [7]:  'dcto':       weather regime active in life cycle within a specified period
                                    after decay of current weather regime (counterpart of onfr)
                [8]:  'dctoID':     id of life cycle of weather regime in dcto
                [9]:  'dctoDATE':   first date of life cycle of weather regime in dcto
                [10]: 'onfr':       weather regime active in life cycle within a specified period
                                    before onset of current weather regime (counterpart of dcto)
                [11]: 'onto':       of secondary importance; counterpart of dcfr
                [12]: 'onfrID':     id of life cycle of weather regime in onfr
                [13]: 'onfromDATE': last date of life cycle of weather regime in onfr
                [14]: 'trfr':       weather regime active in life cycle before current weather regime
                                    becomes active in life cycle (= transition within life cycle vector);
                                    difference to onfr: current weather regime does not need to have
                                    onset at time of trfr but potentially already before (overlapping life cycles); counterpart to trto
                [15]: 'trfrID':     id of life cycle of weather regime in trfr
                [16]: 'trfromDATE': date of transition from trfr to current weather regime
                [17]: 'trto':       weather regime active in life cycle after current weather regime
                                    ends being active in life cycle (= transition within life cycle vector);
                                    difference to dcto: current weather regime does not need to have
                                    decay at time of trto but potentially only afterward (overlapping life cycles); counterpart to trfr
                [18]: 'trtoID':    id of life cycle of weather regime in trto
                [19]: 'trtoDATE':  date of transition from current to trto weather regime
            for no regime, the following characteristic values are given for every
            individual no regime "(pseudo) life cycle" (period of no regime):
                [0]: 'number':     life cycle id (only unique within the corresponding weather regime),
                [1]: 'onset':      date of onset,
                [2]: 'decay':      date of decay,
                [3]: 'duration':   duration of life cycle in hours,
                [4]: 'comes_from': weather regime active in life cycle before current no regime
                [5]: 'ID_from':    id of weather regime in comes_from
                [6]: 'transition_to': weather regime active in life cycle after current no regime
                [7]: 'ID_to':      id of weather regime in transition_to

    Author: Dominik Bueeler

    Date: September 2019

    '''

    ##################
    # Initialization #
    ##################

    # Import
    import sys
    import numpy as np
    import datetime as dt
    from collections import OrderedDict as odict
    import os

      

    # Set basic filenames
    if dataset == 'erainterim':
        if setup == 'z500anom_1979_2015_on_wrdef_10d_1.0_1979_2015': # 10d low-pass filter and iwr threshold of 1.0
            basefname = 'Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local_cluster'
    elif dataset == 'era5':
        if setup == 'z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019': # 10d low-pass filter and iwr threshold of 1.0
            basefname = 'Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local_cluster'

    # Set era-interim reference infile (to get weather regime index and name assignment)
    if dataset == 'erainterim':
        if setup == 'z500anom_1979_2015_on_wrdef_10d_1.0_1979_2015': # 10d low-pass filter and iwr threshold of 1.0
            basefname_iwr_ref = 'Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_proj_local'
            infile_iwr_ref    = '%s/%s.txt' % (basepath, basefname_iwr_ref)
    elif dataset == 'era5':
        if setup == 'z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019': # 10d low-pass filter and iwr threshold of 1.0
            basefname_iwr_ref = 'Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_proj_local'
            infile_iwr_ref    = '%s/%s.txt' % (basepath, basefname_iwr_ref)

    # Get indices of weather regimes from data (do not change!)
    reffile = open(infile_iwr_ref, 'r')
    for l in range(3):
        line = reffile.readline()
    wrsind                       = line.split()[4:]
    if ('ZOEA' in wrsind) and ('ZOWE' in wrsind) and ('BL' in wrsind): #this renames old abbreviations used for erainterim
        wrsind[wrsind.index('ZOEA')] = 'ScTr'
        wrsind[wrsind.index('ZOWE')] = 'EuBL'
        wrsind[wrsind.index('BL')]   = 'ScBL'
   
    # Define weather regime names and old names and corresponding indices (do not change!)
    wrmeta = np.empty(shape=8, dtype=[('wrname','S4'), ('wrname_old','S4'), ('wrindex',int)])
    wrmeta[0]['wrname'], wrmeta[0]['wrname_old'], wrmeta[0]['wrindex'] = 'AT', 'AT', wrsind.index('AT')+1
    wrmeta[1]['wrname'], wrmeta[1]['wrname_old'], wrmeta[1]['wrindex'] = 'ZO', 'ZO', wrsind.index('ZO')+1
    wrmeta[2]['wrname'], wrmeta[2]['wrname_old'], wrmeta[2]['wrindex'] = 'ScTr', 'ZOEA', wrsind.index('ScTr')+1
    wrmeta[3]['wrname'], wrmeta[3]['wrname_old'], wrmeta[3]['wrindex'] = 'AR', 'AR', wrsind.index('AR')+1
    wrmeta[4]['wrname'], wrmeta[4]['wrname_old'], wrmeta[4]['wrindex'] = 'EuBL', 'ZOWE', wrsind.index('EuBL')+1
    wrmeta[5]['wrname'], wrmeta[5]['wrname_old'], wrmeta[5]['wrindex'] = 'ScBL', 'BL', wrsind.index('ScBL')+1
    wrmeta[6]['wrname'], wrmeta[6]['wrname_old'], wrmeta[6]['wrindex'] = 'GL', 'GL', wrsind.index('GL')+1
    wrmeta[7]['wrname'], wrmeta[7]['wrname_old'], wrmeta[7]['wrindex'] = 'no', 'none', 0

    
    #############
    # Load data #
    #############

    # Initialize structured array (temporary)
    data_orig = odict()

    # Loop over weather regimes
    for (w, wr_tmp) in enumerate(wr):

        # Set current weather regime
        wrname     = wrmeta[wrmeta['wrname']==wr_tmp.encode()]['wrname'][0]
        wrname_old = wrmeta[wrmeta['wrname']==wr_tmp.encode()]['wrname_old'][0]
        wrindex    = wrmeta[wrmeta['wrname']==wr_tmp.encode()]['wrindex'][0]

        # Set current infile
        infile = '%s/%s%s.txt' % (basepath, basefname, str(wrindex))

        # Check if weather regime in the file matches with the current weather regime in the code
        reffile = open(infile, 'r')
        for l in range(4):
            line = reffile.readline()
        clusterind_file = int(line.split()[5:][0])
        reffile = open(infile, 'r')
        for l in range(5):
            line = reffile.readline()
        wr_file = line.split()[0]
        # print(infile, wrindex, wrname, clusterind_file, wr_file)
        if (clusterind_file != int(wrindex)) or (dataset=="erainterim" and wr_file != wrname_old.decode()):
            sys.exit('weather regime in the file does not match with the weather regime of the code!')

        # Load life cycle data of current weather regime
        if not wrindex == 0:
 
            if dataset=='erainterim': 
              data_orig[wrname.decode()] = np.loadtxt(infile,
                                                    usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                                                    dtype = [('number','int'),
                                                             ('onset','S11'),
                                                             ('sat_start','S11'),
                                                             ('mx','S11'),
                                                             ('sat_end','S11'),
                                                             ('decay','S11'),
                                                             ('dcfr','S4'),
                                                             ('dcto','S4'),
                                                             ('dctoID','int'),
                                                             ('dctoDATE','S11'),
                                                             ('onfr','S4'),
                                                             ('onto','S4'),
                                                             ('onfrID','int'),
                                                             ('onfromDATE','S11'),
                                                             ('trfr','S4'),
                                                             ('trfrID','int'),
                                                             ('trfromDATE','S11'),
                                                             ('trto','S4'),
                                                             ('trtoID','int'),
                                                             ('trtoDATE','S11')],
                                                    converters = {6: lambda x: wrmeta[wrmeta['wrname_old']==x]['wrname'][0],
                                                                  7: lambda x: wrmeta[wrmeta['wrname_old']==x]['wrname'][0],
                                                                  10: lambda x: wrmeta[wrmeta['wrname_old']==x]['wrname'][0],
                                                                  11: lambda x: wrmeta[wrmeta['wrname_old']==x]['wrname'][0],
                                                                  14: lambda x: wrmeta[wrmeta['wrname_old']==x]['wrname'][0],
                                                                  17: lambda x: wrmeta[wrmeta['wrname_old']==x]['wrname'][0]},
                                                    skiprows = 11)
            elif dataset=='era5': 
              data_orig[wrname.decode()] = np.loadtxt(infile,
                                                    usecols = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
                                                    dtype = [('number','int'),
                                                             ('onset','S11'),
                                                             ('sat_start','S11'),
                                                             ('mx','S11'),
                                                             ('sat_end','S11'),
                                                             ('decay','S11'),
                                                             ('dcfr','S4'),
                                                             ('dcto','S4'),
                                                             ('dctoID','int'),
                                                             ('dctoDATE','S11'),
                                                             ('onfr','S4'),
                                                             ('onto','S4'),
                                                             ('onfrID','int'),
                                                             ('onfromDATE','S11'),
                                                             ('trfr','S4'),
                                                             ('trfrID','int'),
                                                             ('trfromDATE','S11'),
                                                             ('trto','S4'),
                                                             ('trtoID','int'),
                                                             ('trtoDATE','S11')],
                                                    skiprows = 11)
        elif wrindex == 0:
          if dataset=='erainterim':
              data_orig[wrname.decode()] = np.loadtxt(infile,
                                                    usecols = (0,1,2,3,4,5,6,7),
                                                    dtype = [('number','int'),
                                                             ('onset','S11'),
                                                             ('decay','S11'),
                                                             ('duration','int'),
                                                             ('comes_from','S4'),
                                                             ('ID_from','int'),
                                                             ('transition_to','S4'),
                                                             ('ID_to','int')],
                                                    converters = {4: lambda x: wrmeta[wrmeta['wrname_old']==x]['wrname'][0],
                                                                  6: lambda x: wrmeta[wrmeta['wrname_old']==x]['wrname'][0]},
                                                    skiprows = 11)
          elif dataset=='era5':
              data_orig[wrname.decode()] = np.loadtxt(infile,
                                                    usecols = (0,1,2,3,4,5,6,7),
                                                    dtype = [('number','int'),
                                                             ('onset','S11'),
                                                             ('decay','S11'),
                                                             ('duration','int'),
                                                             ('comes_from','S4'),
                                                             ('ID_from','int'),
                                                             ('transition_to','S4'),
                                                             ('ID_to','int')],
                                                    skiprows = 11)
    
    ####################
    # Postprocess data #
    ####################

    # If time strings have to be converted to datetime objects
    if tformat == 'dtime':

        # Initialize structured array
        data = odict()

        # Loop over weather regimes
        for (w, wr_tmp) in enumerate(wr):

            # Set current weather regime
            wrname  = wrmeta[wrmeta['wrname']==wr_tmp.encode()]['wrname'][0]
            wrindex = wrmeta[wrmeta['wrname']==wr_tmp.encode()]['wrindex'][0]

            # If weather regime is not no regime
            if not wrindex == 0:

                # Initialize structured array for current weather regime
                data[wrname.decode()] = np.empty(shape = data_orig[wrname.decode()].shape,
                                                 dtype = [('number','int'),
                                                          ('onset',object),
                                                          ('sat_start',object),
                                                          ('mx',object),
                                                          ('sat_end',object),
                                                          ('decay',object),
                                                          ('dcfr','S4'),
                                                          ('dcto','S4'),
                                                          ('dctoID','int'),
                                                          ('dctoDATE',object),
                                                          ('onfr','S4'),
                                                          ('onto','S4'),
                                                          ('onfrID','int'),
                                                          ('onfromDATE',object),
                                                          ('trfr','S4'),
                                                          ('trfrID','int'),
                                                          ('trfromDATE',object),
                                                          ('trto','S4'),
                                                          ('trtoID','int'),
                                                          ('trtoDATE',object)])

                # Append data to structured array with new date format
                data[wrname.decode()]['number']     = data_orig[wrname.decode()]['number']
                data[wrname.decode()]['onset']      = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['onset'] ])
                data[wrname.decode()]['sat_start']  = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['sat_start'] ])
                data[wrname.decode()]['mx']         = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['mx'] ])
                data[wrname.decode()]['sat_end']    = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['sat_end'] ])
                data[wrname.decode()]['decay']      = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['decay'] ])
                data[wrname.decode()]['dcfr']       = data_orig[wrname.decode()]['dcfr']
                data[wrname.decode()]['dcto']       = data_orig[wrname.decode()]['dcto']
                data[wrname.decode()]['dctoID']     = data_orig[wrname.decode()]['dctoID']
                data[wrname.decode()]['dctoDATE']   = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['dctoDATE'] ])
                data[wrname.decode()]['onfr']       = data_orig[wrname.decode()]['onfr']
                data[wrname.decode()]['onto']       = data_orig[wrname.decode()]['onto']
                data[wrname.decode()]['onfrID']     = data_orig[wrname.decode()]['onfrID']
                data[wrname.decode()]['onfromDATE'] = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['onfromDATE'] ])
                data[wrname.decode()]['trfr']       = data_orig[wrname.decode()]['trfr']
                data[wrname.decode()]['trfrID']     = data_orig[wrname.decode()]['trfrID']
                data[wrname.decode()]['trfromDATE'] = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['trfromDATE'] ])
                data[wrname.decode()]['trto']       = data_orig[wrname.decode()]['trto']
                data[wrname.decode()]['trtoID']     = data_orig[wrname.decode()]['trtoID']
                data[wrname.decode()]['trtoDATE']   = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                 for t in data_orig[wrname.decode()]['trtoDATE'] ])

            # If weather regime is no regime
            elif wrindex == 0:

                # Initialize structured array for current weather regime
                data[wrname.decode()] = np.empty(shape = data_orig[wrname.decode()].shape,
                                                 dtype = [('number','int'),
                                                          ('onset',object),
                                                          ('decay',object),
                                                          ('duration','int'),
                                                          ('comes_from','S4'),
                                                          ('ID_from','int'),
                                                          ('transition_to','S4'),
                                                          ('ID_to','int')])

                # Append data to structured array
                data[wrname.decode()]['number']        = data_orig[wrname.decode()]['number']
                data[wrname.decode()]['onset']         = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                    for t in data_orig[wrname.decode()]['onset'] ])
                data[wrname.decode()]['decay']         = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                                                    for t in data_orig[wrname.decode()]['decay'] ])
                data[wrname.decode()]['duration']      = data_orig[wrname.decode()]['duration']
                data[wrname.decode()]['comes_from']    = data_orig[wrname.decode()]['comes_from']
                data[wrname.decode()]['ID_from']       = data_orig[wrname.decode()]['ID_from']
                data[wrname.decode()]['transition_to'] = data_orig[wrname.decode()]['transition_to']
                data[wrname.decode()]['ID_to']         = data_orig[wrname.decode()]['ID_to']

    # Loop over weather regimes
    for (w, wr_tmp) in enumerate(wr):

        # Set current weather regime
        wrname = wrmeta[wrmeta['wrname']==wr_tmp.encode()]['wrname'][0]

        # If time should be returned as strings
        if tformat == 'string':

            # Find index for first onset within specified period
            times_onset = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                     for t in data_orig[wrname.decode()]['onset'] ])
            diff_start  = times_onset - dt.datetime.strptime(start,'%Y%m%d_%H')
            diff_start  = np.array([ d.days for d in diff_start ])
            if np.any(diff_start < 0):
                ind_start = np.where(diff_start < 0)[0][-1] + 1
            else:
                ind_start = 0
                
            # Find index for last decay within specified period
            times_decay = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H')
                                     for t in data_orig[wrname.decode()]['decay'] ])
            diff_end    = times_decay - dt.datetime.strptime(end,'%Y%m%d_%H')
            diff_end    = np.array([ d.days for d in diff_end ])
            if np.any(diff_end > 0):
                ind_end = np.where(diff_end > 0)[0][0] - 1
            else:
                ind_end = len(data_orig[wrname.decode()])

            # Trim data according to specified period
            data_orig[wrname.decode()] = data_orig[wrname.decode()][ind_start:ind_end+1]

        # If time should be returned as datetime objects
        elif tformat == 'dtime':

            # Find index for first onset within specified period
            times_onset = data[wrname.decode()]['onset']
            diff_start  = times_onset - dt.datetime.strptime(start,'%Y%m%d_%H')
            diff_start  = np.array([ d.days for d in diff_start ])
            if np.any(diff_start < 0):
                ind_start = np.where(diff_start < 0)[0][-1] + 1
            else:
                ind_start = 0

            # Find index for last decay within specified period
            times_decay = data[wrname.decode()]['decay']
            diff_end    = times_decay - dt.datetime.strptime(end,'%Y%m%d_%H')
            diff_end    = np.array([ d.days for d in diff_end ])
            if np.any(diff_end > 0):
                ind_end = np.where(diff_end > 0)[0][0] - 1
            else:
                ind_end = len(data[wrname.decode()])
                
            # Trim data according to specified period
            data[wrname.decode()] = data[wrname.decode()][ind_start:ind_end+1]

    # Return
    if tformat == 'string':
        return data_orig
    elif tformat == 'dtime':
        return data

####################
# Execute function #
####################

# If function is directly executed
if __name__ == '__main__':

    # Import
    from fct_wrlcera import wrlcera

    # Execute function
    data = wrlcera()

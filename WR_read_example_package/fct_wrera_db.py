#!/usr/bin/env python

# original code by Dominik Büelerdominik.bueeler@meteoswiss.ch
# 20250102: Simplfications by christian.grams@gmx.de
#          


###################
# Define function #
###################

# Function
def wrera(start, end, hours, tformat, setup, dataset, basepath):

    '''
    Description: Load weather regime indices and life cycles from era-interim or era5
                 for a specific period (for details about the weather regimes
                 contact Christian Grams, christian.grams@gmx.de)

    Input: start   = start date as string in the format YYYYMMDD_HH
           end     = end date as string in the format YYYYMMDD_HH
           hours   = hour or list of hours ('00','03','06','09','12','15','18','21') which should be included
           tformat = 'string' if dates should be returned in string format,
                     'dtime' if dates should be returned in datetime format
           setup = string defining the weather regime setup from which the data should
                   be read -> it is important to select this carefully!;
                   'z500anom_1979_2015_on_wrdef_10d_1.0_1979_2015' = z500 anomalies with respect
                   to climatology of 1979-2015 are projected on weather regime cluster means
                   based on 10d-low-pass-filtered z500 anomalies between 1979 and 2015
                   with an iwr threshold of 1.0 (erainterim),
                   'z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019' = z500 anomalies with respect
                   to climatology of 1979-2019 are projected on weather regime cluster means
                   based on 10d-low-pass-filtered z500 anomalies between 1979 and 2019
                   with an iwr threshold of 1.0 (era5, default setup),
           dataset  = 'erainterim' or 'era5'
           basepath = directory to the raw data
           
    Output: Dictionary with the general structure dict[key], with the keys being the following:
            'IWR': weather regime index after michel & riviere (2011) with the fields
                   time since 19790101_00 in hours ('tsince'), valid time ('time'),
                   cluster class index ('cci'), and the different
                   weather regime projections ('AT','ZO','ScTr','AR','EuBL','ScBL','GL');
                   note that order of weather regimes in files can change depending on dataset!
            'MAXIWR': maximum weather regime index after michel & riviere (2011)
                      with the fields time since 19790101_00 in hours ('tsince'), valid time ('time'),
                      index of maximum weather regime projection ('wrindex'),
                      and name of maximum weather regime projection ('wrname');
                      note that assignment of wrindex to wrname can change depending on dataset!
            'LC': full life cycle with the fields time since ('tsince') 19790101_00 in hours,
                  valid time ('time'), index of life cycle ('wrindex'), and name of life cycle ('wrname');
                  note that assignment of wrindex to wrname can change depending on dataset!

    Author: Dominik Bueeler

    Date: November 2020

    '''

    ##################
    # Initialization #
    ##################

    # Import
    import numpy as np
    import datetime as dt
    from collections import OrderedDict as odict
    import os
    

    
    # Set basic filenames
    if dataset == 'erainterim':
        if setup == 'z500anom_1979_2015_on_wrdef_10d_1.0_1979_2015': # 10d low-pass filter and iwr threshold of 1.0
            basefname_iwr    = 'Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_proj_local'
            basefname_maxiwr = 'Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local'
            basefname_lc     = 'Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local'
    elif dataset == 'era5':
        if setup == 'z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019': # 10d low-pass filter and iwr threshold of 1.0
            basefname_iwr    = 'Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_proj_local'
            basefname_maxiwr = 'Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local'
            basefname_lc     = 'Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local'

    # Set era-interim infile
    infile_iwr    = '%s/%s.txt' % (basepath, basefname_iwr)
    infile_maxiwr = '%s/%s.txt' % (basepath, basefname_maxiwr)
    infile_lc     = '%s/%s.txt' % (basepath, basefname_lc)

    # Get order of weather regimes from data (do not change!)
    reffile = open(infile_iwr, 'r')
    for l in range(5):
        line = reffile.readline()
    wrsorder = line.split()[16:]
    if ('ZOEA' in wrsorder) and ('ZOWE' in wrsorder) and ('BL' in wrsorder):
        wrsorder[wrsorder.index('ZOEA')] = 'ScTr'
        wrsorder[wrsorder.index('ZOWE')] = 'EuBL'
        wrsorder[wrsorder.index('BL')]   = 'ScBL'

    # Get indices of weather regimes from data (do not change!)
    reffile = open(infile_iwr, 'r')
    for l in range(3):
        line = reffile.readline()
    wrsind = line.split()[4:]
    if ('ZOEA' in wrsind) and ('ZOWE' in wrsind) and ('BL' in wrsind):
        wrsind[wrsind.index('ZOEA')] = 'ScTr'
        wrsind[wrsind.index('ZOWE')] = 'EuBL'
        wrsind[wrsind.index('BL')]   = 'ScBL'

    # Define weather regime indices and corresponding names (do not change!)
    wrmeta = np.empty(shape=8, dtype=[('wrname','S4'), ('wrindex',int)])
    wrmeta[0]['wrname'], wrmeta[0]['wrindex'] = 'AT', wrsind.index('AT')+1
    wrmeta[1]['wrname'], wrmeta[1]['wrindex'] = 'ZO',  wrsind.index('ZO')+1
    wrmeta[2]['wrname'], wrmeta[2]['wrindex'] = 'ScTr', wrsind.index('ScTr')+1
    wrmeta[3]['wrname'], wrmeta[3]['wrindex'] = 'AR', wrsind.index('AR')+1
    wrmeta[4]['wrname'], wrmeta[4]['wrindex'] = 'EuBL', wrsind.index('EuBL')+1
    wrmeta[5]['wrname'], wrmeta[5]['wrindex'] = 'ScBL', wrsind.index('ScBL')+1
    wrmeta[6]['wrname'], wrmeta[6]['wrindex'] = 'GL', wrsind.index('GL')+1
    wrmeta[7]['wrname'], wrmeta[7]['wrindex'] = 'no', 0

    #############
    # Load data #
    #############

    # Initialize data types
    dtype_iwr              = [('tsince','int'),
                              ('time','S11'),
                              ('cci','int')] + list(zip(wrsorder,['float']*7))
    dtype_iwr_dtime        = [('tsince','int'),
                              ('time',object),
                              ('cci','int')] + list(zip(wrsorder,['float']*7))
    dtype_maxiwr           = [('tsince','int'),
                              ('time','S11'),
                              ('wrindex','int')]
    dtype_maxiwr_ext       = [('tsince','int'),
                              ('time','S11'),
                              ('wrindex','int'),
                              ('wrname','S4')]
    dtype_maxiwr_ext_dtime = [('tsince','int'),
                              ('time',object),
                              ('wrindex','int'),
                              ('wrname','S4')]
    dtype_lc               = [('tsince','int'),
                              ('time','S11'),
                              ('wrindex','int')]
    dtype_lc_ext           = [('tsince','int'),
                              ('time','S11'),
                              ('wrindex','int'),
                              ('wrname','S4')]
    dtype_lc_ext_dtime     = [('tsince','int'),
                              ('time',object),
                              ('wrindex','int'),
                              ('wrname','S4')]

    # Initialize original structured array (temporary)
    data_orig = odict()

    # Load era-interim data
    data_orig['IWR']    = np.loadtxt(infile_iwr,
                                     usecols = (0,1,2,4,5,6,7,8,9,10),
                                     dtype = dtype_iwr,
                                     skiprows = 7)
    data_orig['MAXIWR'] = np.loadtxt(infile_maxiwr,
                                     usecols = (0,1,3),
                                     dtype = dtype_maxiwr,
                                     skiprows = 7)
    data_orig['LC']     = np.loadtxt(infile_lc,
                                     usecols = (0,1,4),
                                     dtype = dtype_lc,
                                     skiprows = 7)
    
    ####################
    # Postprocess data #
    ####################

    # Convert datestring to datetime vector
    dtimes = np.array([ dt.datetime.strptime(t.decode(),'%Y%m%d_%H') for t in data_orig['IWR']['time'] ])
    
    # If time strings do not have to be converted to datetime objects
    if tformat == 'string':

        # Initialize output structured arrays
        data           = odict()
        data['IWR']    = np.empty(shape = data_orig['IWR'].shape,
                                  dtype = dtype_iwr)
        data['MAXIWR'] = np.empty(shape = data_orig['MAXIWR'].shape,
                                  dtype = dtype_maxiwr_ext)
        data['LC']     = np.empty(shape = data_orig['LC'].shape,
                                  dtype = dtype_lc_ext)

        # Append data to output structured arrays
        data['IWR']['tsince']     = data_orig['IWR']['tsince']
        data['IWR']['time']       = data_orig['IWR']['time']
        data['IWR']['cci']        = data_orig['IWR']['cci']
        data['IWR']['AT']         = data_orig['IWR']['AT']
        data['IWR']['ZO']         = data_orig['IWR']['ZO']
        data['IWR']['ScTr']       = data_orig['IWR']['ScTr']
        data['IWR']['AR']         = data_orig['IWR']['AR']
        data['IWR']['EuBL']       = data_orig['IWR']['EuBL']
        data['IWR']['ScBL']       = data_orig['IWR']['ScBL']
        data['IWR']['GL']         = data_orig['IWR']['GL']
        data['MAXIWR']['tsince']  = data_orig['MAXIWR']['tsince']
        data['MAXIWR']['time']    = data_orig['MAXIWR']['time']
        data['MAXIWR']['wrindex'] = data_orig['MAXIWR']['wrindex']
        data['MAXIWR']['wrname']  = np.array([ wrmeta[wrmeta['wrindex']==i]['wrname'][0]
                                               for i in data_orig['MAXIWR']['wrindex'] ])
        data['LC']['tsince']      = data_orig['LC']['tsince']
        data['LC']['time']        = data_orig['LC']['time']
        data['LC']['wrindex']     = data_orig['LC']['wrindex']
        data['LC']['wrname']      = np.array([ wrmeta[wrmeta['wrindex']==i]['wrname'][0]
                                               for i in data_orig['LC']['wrindex'] ])

    # If time strings have to be converted to datetime objects
    elif tformat == 'dtime':

        # Initialize output structured arrays
        data           = odict()
        data['IWR']    = np.empty(shape = data_orig['IWR'].shape,
                                  dtype = dtype_iwr_dtime)
        data['MAXIWR'] = np.empty(shape = data_orig['MAXIWR'].shape,
                                  dtype = dtype_maxiwr_ext_dtime)
        data['LC']     = np.empty(shape = data_orig['LC'].shape,
                                  dtype = dtype_lc_ext_dtime)

        # Append data to output structured arrays with new date format
        data['IWR']['tsince']     = data_orig['IWR']['tsince']
        data['IWR']['time']       = dtimes
        data['IWR']['cci']        = data_orig['IWR']['cci']
        data['IWR']['AT']         = data_orig['IWR']['AT']
        data['IWR']['ZO']         = data_orig['IWR']['ZO']
        data['IWR']['ScTr']       = data_orig['IWR']['ScTr']
        data['IWR']['AR']         = data_orig['IWR']['AR']
        data['IWR']['EuBL']       = data_orig['IWR']['EuBL']
        data['IWR']['ScBL']       = data_orig['IWR']['ScBL']
        data['IWR']['GL']         = data_orig['IWR']['GL']
        data['MAXIWR']['tsince']  = data_orig['MAXIWR']['tsince']
        data['MAXIWR']['time']    = dtimes
        data['MAXIWR']['wrindex'] = data_orig['MAXIWR']['wrindex']
        data['MAXIWR']['wrname']  = np.array([ wrmeta[wrmeta['wrindex']==i]['wrname'][0]
                                               for i in data_orig['MAXIWR']['wrindex'] ])
        data['LC']['tsince']      = data_orig['LC']['tsince']
        data['LC']['time']        = dtimes
        data['LC']['wrindex']     = data_orig['LC']['wrindex']
        data['LC']['wrname']      = np.array([ wrmeta[wrmeta['wrindex']==i]['wrname'][0]
                                               for i in data_orig['LC']['wrindex'] ])

    # Get indices of desired start and end times
    if tformat == 'string':
        ind_start = np.where(data['IWR']['time'] == start.encode())[0][0]
        ind_end   = np.where(data['IWR']['time'] == end.encode())[0][0]
    elif tformat == 'dtime':
        ind_start = np.where(data['IWR']['time'] == dt.datetime.strptime(start,'%Y%m%d_%H'))[0][0]
        ind_end   = np.where(data['IWR']['time'] == dt.datetime.strptime(end,'%Y%m%d_%H'))[0][0]

    # Trim data according to time period
    data['IWR']    = data['IWR'][ind_start:ind_end+1]
    data['MAXIWR'] = data['MAXIWR'][ind_start:ind_end+1]
    data['LC']     = data['LC'][ind_start:ind_end+1]
    dtimes         = dtimes[ind_start:ind_end+1]

    # Get indices of desired hours
    if dataset == 'erainterim':
        if not hours == ['00','06','12','18']:
            ind_00 = []
            ind_06 = []
            ind_12 = []
            ind_18 = []
            if hours == '00' or '00' in hours:
                dtimes_00 = np.array([ t for t in dtimes if t.hour == 0 ])
                ind_00    = [ np.where(dtimes == t)[0][0] for t in dtimes_00 ]
            if hours == '06' or '06' in hours:
                dtimes_06 = np.array([ t for t in dtimes if t.hour == 6 ])
                ind_06    = [ np.where(dtimes == t)[0][0] for t in dtimes_06 ]
            if hours == '12' or '12' in hours:
                dtimes_12 = np.array([ t for t in dtimes if t.hour == 12 ])
                ind_12    = [ np.where(dtimes == t)[0][0] for t in dtimes_12 ]
            if hours == '18' or '18' in hours:
                dtimes_18 = np.array([ t for t in dtimes if t.hour == 18 ])
                ind_18    = [ np.where(dtimes == t)[0][0] for t in dtimes_18 ]
            ind = sorted(ind_00 + ind_06 + ind_12 + ind_18)
    elif dataset == 'era5':
        if not hours == ['00','03','06','09','12','15','18','21']:
            ind_00 = []
            ind_03 = []
            ind_06 = []
            ind_09 = []
            ind_12 = []
            ind_15 = []
            ind_18 = []
            ind_21 = []
            if hours == '00' or '00' in hours:
                dtimes_00 = np.array([ t for t in dtimes if t.hour == 0 ])
                ind_00    = [ np.where(dtimes == t)[0][0] for t in dtimes_00 ]
            if hours == '03' or '03' in hours:
                dtimes_03 = np.array([ t for t in dtimes if t.hour == 3 ])
                ind_03    = [ np.where(dtimes == t)[0][0] for t in dtimes_03 ]
            if hours == '06' or '06' in hours:
                dtimes_06 = np.array([ t for t in dtimes if t.hour == 6 ])
                ind_06    = [ np.where(dtimes == t)[0][0] for t in dtimes_06 ]
            if hours == '09' or '09' in hours:
                dtimes_09 = np.array([ t for t in dtimes if t.hour == 9 ])
                ind_09    = [ np.where(dtimes == t)[0][0] for t in dtimes_09 ]
            if hours == '12' or '12' in hours:
                dtimes_12 = np.array([ t for t in dtimes if t.hour == 12 ])
                ind_12    = [ np.where(dtimes == t)[0][0] for t in dtimes_12 ]
            if hours == '15' or '15' in hours:
                dtimes_15 = np.array([ t for t in dtimes if t.hour == 15 ])
                ind_15    = [ np.where(dtimes == t)[0][0] for t in dtimes_15 ]
            if hours == '18' or '18' in hours:
                dtimes_18 = np.array([ t for t in dtimes if t.hour == 18 ])
                ind_18    = [ np.where(dtimes == t)[0][0] for t in dtimes_18 ]
            if hours == '21' or '21' in hours:
                dtimes_21 = np.array([ t for t in dtimes if t.hour == 21 ])
                ind_21    = [ np.where(dtimes == t)[0][0] for t in dtimes_21 ]
            ind = sorted(ind_00 + ind_03 + ind_06 + ind_09 + ind_12 + ind_15 + ind_18 + ind_21)

    # Extract data with desired hours
    if dataset == 'erainterim':
        if not hours == ['00','06','12','18']:
            data['LC']     = data['LC'][ind]
            data['MAXIWR'] = data['MAXIWR'][ind]
            data['IWR']    = data['IWR'][ind]
    elif dataset == 'era5':
        if not hours == ['00','03','06','09','12','15','18','21']:
            data['LC']     = data['LC'][ind]
            data['MAXIWR'] = data['MAXIWR'][ind]
            data['IWR']    = data['IWR'][ind]
        
    # Return
    return dtimes, data

####################
# Execute function #
####################

# If function is directly executed
if __name__ == '__main__':

    # Import
    from fct_wrera import wrera

    # Define input parameters
    # start   = '19790111_00'
    # end     = '20190121_18'                                                               
    # hours   = ['00','06','12','18']                                 
    # tformat = 'string'                                                                
    # setup   = 'z500anom_1979_2018_on_wrdef_5d_0.9_1979_2018'                           
    # dataset = 'erainterim'                                                                  
    # env     = 'local'                                                                     
    start   = '19790111_00'
    end     = '20191231_21'                                                               
    hours   = ['00','03','06','09','12','15','18','21']                                 
    tformat = 'string'                                                                
    setup   = 'z500anom_1979_2019_on_wrdef_10d_1.0_1979_2019'                           
    dataset = 'era5'                                                                  
    basepath = './ec.era5/eof/update_1950/latwgt/eof_AUTO_Z500_N161_-80E40E30N90N_year/norminput'                                                                     

    # Execute function
    data = wrera(start = start,
                 end = end,
                 hours = hours,
                 tformat = tformat,
                 setup = setup,
                 dataset = dataset,
                 basepath = basepath)

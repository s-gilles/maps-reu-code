This is the code that was used to generate the results of the 2014
MAPS REU (See http://www.curve.unhyperbolic.org/activities.html and
http://www-math.umd.edu/maps-reu.html ). In general:

 - `all_volumes.csv` is a compilation of the result of

        #!/usr/bin/python2 -u
        from snappy import *
        from VolumeFinder import *
    
        li = iter(OrientableCuspedCensus[1234:5678])
    
        # Iterators customized to taste
        begin_collection(BatchIterator(DehnFillIterator(li), 2500),
                                       output_filename = 'partial.csv',
                                       thread_num = 12,
                                       is_appending = False)
    

 - `spans.csv` is a compilation of the results of

        #!/usr/bin/python2 -u
        from VolumeProcessing import *
    
        quick_write_spans([
            'partial_dehn_fills_00000-17000.csv',
            'partial_orientablecuspedcensus_1.csv',
            # ... ,
            'some_other_files.csv'],
            'output_volume_spans.csv',
            borel_shape_field_degree = 48, skip_borel = False)
    

 - `linear_combinations_...csv` are results of running
   `linear_combinations_master.py` with appropriate options (this one
   respects `--help`).

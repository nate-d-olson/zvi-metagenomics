Scripts to perform sanity checks on input data (e.g. OTU tables), sample metadata, and primary anlysis output. 

__Metadata__ 
checks: correct number of samples and values for all fields 

__Input data__ 
* OTU table
    - checks: metadata sample names and column names match, no negative values  

__Primary analysis__
* data frame with alpha diversity metrics    
    - checks: diversity metric values for all samples, no negative values  
* matrix for each beta diversity metric  
    - checks: number of columns and row equal to sample number, all values between 0 and 1 

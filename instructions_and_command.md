### Access Bridges-2 Server.
'''
    $ ssh jji5@bridges2.psc.edu
'''
Then enter password for your bridges-2 account. You won't be able to see the password when you type. 

### Find data storing folder.
'''
    $ cd /ocean/projects/bio230007p/ikaplow/HumanAtac/Pancreas/peak
'''
There are many subfolders inside. For now we should use idr_reproducibility
'''
    $ cd idr_reproducibility
'''
To look at the data structure of one of the files, we can do
'''
    $ zcat idr.optimal_peak.narrowPeak.gz | head
'''
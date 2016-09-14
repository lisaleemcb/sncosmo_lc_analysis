#!/usr/bin/anv python

import analyzeSN as ans


def test_loads_full():
    """
    Test that the entire fits file can be loaded correctly. We only do
    a partial check on 'correctness' by checking that
    1. the number of SN is correct
    2. The SNID is correct
    3. The number of epochs is correct
    by reading off the number in the example file.
    """
    sne = ans.snanaSims.SnanaSims.fromSNANAfileroot('snana_fits',
                                                    location=ans.example_data)
    snlist = sne.snList
    # Test 1.
    assert len(snlist) == 2
    # Test 2.
    assert snlist[0].meta['SNID'] == '03D1aw'
    assert snlist[1].meta['SNID'] == '03D1ax'
    # Test 3.
    assert len(snlist[0]) == 48


def test_loads_bysnid():
    """
    Test that the fits file can be loaded correctly by a list of SNIDs.
    We only do
    a partial check on 'correctness' by checking that
    1. the number of SN is correct
    2. The SNID is correct
    3. The number of epochs is correct
    by reading off the number in the example file.
    """

    sne = ans.snanaSims.SnanaSims.fromSNANAfileroot('snana_fits',
                                                    location=ans.example_data,
                                                    snids=['03D1aw'])
    snlist = sne.snList
    assert len(snlist) == 1
    assert snlist[0].meta['SNID'] == '03D1aw'
    assert len(snlist[0]) == 48
    
def test_loads_bynum():
    """
    Test that the fits file can be loaded correctly by the first n objects
    We only do
    a partial check on 'correctness' by checking that
    1. the number of SN is correct
    2. The SNID is correct
    3. The number of epochs is correct
    by reading off the number in the example file.
    """

    sne = ans.snanaSims.SnanaSims.fromSNANAfileroot('snana_fits',
                                                    location=ans.example_data,
                                                    n=1)
    snlist = sne.snList
    assert len(snlist) == 1
    assert snlist[0].meta['SNID'] == '03D1aw'
    assert len(snlist[0]) == 48
if __name__ == '__main__':
    test_loads_full()

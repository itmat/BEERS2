import functools

def demultiplexer(lookup_dict, max_allowed_errors=1):
    '''
    Demultiplexer based off a barcode allowing some number of mismatch errors
    Finds the appropriate object in lookup_dict (which maps exact barcodes to some objects,
    and should be a default dict). If demultiplexing fails, then use the default object of the dict.

    Returns a function which performs the demultiplexing, i.e. call this as
    demux = demultiplexer({"ACG": thing1, "CCG": thing2})
    demux("ACT") # yields thing1
    '''
    # Grab the keys now since they will change as things are looked up
    # default dict adds keys when they are looked up
    keys = list(lookup_dict.keys())

    # We make a cache for efficiency since this lookup is done on every single molecule
    # and since lookup_dict can't be cached hashed for the cache, we have to make a closure
    # here wrapping around it
    @functools.lru_cache(maxsize=None)
    def demultiplex(barcode):
        def num_mismatches(barcode, test):
            return sum(1 if b != t else 0 for b,t in zip(barcode, test))
        matches = [k for k in keys if num_mismatches(barcode, k) <= max_allowed_errors]
        assert len(matches) <= 1, f"Too many matches {matches} for {barcode}" # Barcodes should be chosen so that none are within max_allowed_errors of each other
        if len(matches) == 0:
            return lookup_dict[barcode] # Will use the default object
        else:
            return lookup_dict[matches[0]] # Use the match (possibly exact)
    return demultiplex


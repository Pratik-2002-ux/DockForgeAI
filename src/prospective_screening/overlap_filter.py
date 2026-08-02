def exclude_overlap(frame,training_keys,external_keys): return frame[~frame.inchikey_connectivity.isin(set(training_keys)|set(external_keys))].copy()

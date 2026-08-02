"""Source table merger."""
import pandas as pd
def merge_measurements(*frames): return pd.concat([f for f in frames if f is not None and not f.empty],ignore_index=True)

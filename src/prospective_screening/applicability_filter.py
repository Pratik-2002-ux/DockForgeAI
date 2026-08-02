def in_domain(frame): return frame[frame["applicability_domain"].eq("in_domain")].copy()

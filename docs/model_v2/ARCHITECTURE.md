# DockForgeAI v2 architecture

v2 is isolated under `src/models_v2`, `scripts/model_v2`, `app/components/model_v2`, `data/model_v2`, `reports/model_v2`, and `docs/model_v2`. It does not import v1 trained artifacts or modify v1 outputs.

The dataset schema permits only experimentally measured active and inactive records from ChEMBL, PubChem BioAssay, or BindingDB exports. DUD-E records are rejected by construction. Every experiment records seed, package versions, dataset checksum, Git hash, hyperparameters, runtime, and hardware.

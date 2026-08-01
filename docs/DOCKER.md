# Docker

Start the local application with:

```bash
docker compose up --build
```

Open <http://localhost:8501>. Compose maps `data`, `reports`, and `logs` into the non-root container. The image uses micromamba on Linux, exposes only port 8501, and includes a Streamlit healthcheck.

Apple Silicon runs the Linux/arm64 variant when all conda packages are available. Vina, RDKit, Meeko, Open Babel, and Gemmi availability can differ across architectures. `docker compose config` validates configuration but is not evidence that every image layer built or that docking was executed. Pin or select `platform: linux/amd64` only after accepting emulation cost on constrained hardware.

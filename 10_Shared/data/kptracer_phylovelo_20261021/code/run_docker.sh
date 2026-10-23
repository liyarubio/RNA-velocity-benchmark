#!/usr/bin/env bash
set -euo pipefail
TASK_DIR="$(cd -- "$(dirname -- "$0")/.." && pwd)"
INPUT_DIR="$(dirname -- "$TASK_DIR")/velogif_XC/real_data"
CONTAINER_NAME=velogif-kptracer-phylovelo-20261021
test -s "$TASK_DIR/inputs/3726_NT_T1_tree.nwk"
test -s "$TASK_DIR/inputs/adata_processed.nt.h5ad"
test -s "$TASK_DIR/records/download.json"
if ! docker container inspect "$CONTAINER_NAME" >/dev/null 2>&1; then
 docker run -d --name "$CONTAINER_NAME" --cpus 4 --memory 16g --network none \
  --cap-drop ALL --security-opt no-new-privileges --user "$(id -u):$(id -g)" \
  -e OMP_NUM_THREADS=2 -e OPENBLAS_NUM_THREADS=2 -e MPLCONFIGDIR=/tmp/phylovelo-mpl \
  -e NUMBA_CACHE_DIR=/tmp/phylovelo-numba \
  -v "$TASK_DIR:/task" -v "$INPUT_DIR:/input:ro" \
  --entrypoint /bin/sh velogif:v1 -c 'exec tail -f /dev/null'
else
 docker start "$CONTAINER_NAME"
fi
docker inspect "$CONTAINER_NAME" --format '{{.Image}}' > "$TASK_DIR/records/docker_image_id.txt"
docker exec "$CONTAINER_NAME" /opt/conda/envs/phylovelo/bin/python -u /task/code/run_phylovelo.py > "$TASK_DIR/logs/run.log" 2>&1
docker exec "$CONTAINER_NAME" /opt/conda/envs/phylovelo/bin/python -u /task/code/plot_results.py > "$TASK_DIR/logs/plot.log" 2>&1
echo 'Docker inference and plotting completed'

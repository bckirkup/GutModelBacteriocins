#!/usr/bin/env bash
# LOCAL ONLY — build CUDA image and push to ECR (us-east-1).
# CloudShell cannot build this image. Run from your laptop at the repo root:
#
#   bash deploy/aws/01_push_image.sh
#
# Needs: Docker running, AWS CLI logged in, internet.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/env.sh" ]]; then
  echo "ERROR: env.sh is missing beside this script. Clone the repository and run it from the clone." >&2
  exit 1
fi
# shellcheck source=env.sh
source "${SCRIPT_DIR}/env.sh"

if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: docker not found. Install Docker Desktop and start it." >&2
  exit 1
fi

cd "${ROOT}"

GUTIBM_GIT_SHA="$(git rev-parse HEAD 2>/dev/null || true)"
if [[ -z "${GUTIBM_GIT_SHA}" ]]; then
  GUTIBM_GIT_SHA="unknown"
fi

echo "==> Ensure ECR repository exists"
aws ecr describe-repositories --repository-names gutibm --region "${AWS_REGION}" >/dev/null 2>&1 \
  || aws ecr create-repository --repository-name gutibm --region "${AWS_REGION}" >/dev/null

echo "==> Docker login to ECR"
aws ecr get-login-password --region "${AWS_REGION}" \
  | docker login --username AWS --password-stdin "${ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com"

echo "==> docker build (this can take a long time the first time)"
docker build -f deploy/aws/Dockerfile -t gutibm:cuda \
  --build-arg 'CUDA_ARCHS=75;86;89' \
  --build-arg "GUTIBM_GIT_SHA=${GUTIBM_GIT_SHA}" .

echo "==> tag + push ${IMAGE_URI}"
docker tag gutibm:cuda "${IMAGE_URI}"
docker push "${IMAGE_URI}"

echo "OK: pushed ${IMAGE_URI}"

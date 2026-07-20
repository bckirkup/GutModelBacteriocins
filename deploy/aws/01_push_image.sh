#!/usr/bin/env bash
# LOCAL ONLY — build CUDA image and push to ECR (us-east-1).
# CloudShell cannot build this image. Run from your laptop at the repo root:
#
#   bash deploy/aws/01_push_image.sh
#
# Needs: Docker running, AWS CLI logged in, internet.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=env.sh
source "${ROOT}/deploy/aws/env.sh"

if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: docker not found. Install Docker Desktop and start it." >&2
  exit 1
fi

cd "${ROOT}"

echo "==> Ensure ECR repository exists"
aws ecr describe-repositories --repository-names gutibm --region "${AWS_REGION}" >/dev/null 2>&1 \
  || aws ecr create-repository --repository-name gutibm --region "${AWS_REGION}" >/dev/null

echo "==> Docker login to ECR"
aws ecr get-login-password --region "${AWS_REGION}" \
  | docker login --username AWS --password-stdin "${ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com"

echo "==> docker build (this can take a long time the first time)"
docker build -f deploy/aws/Dockerfile -t gutibm:cuda --build-arg 'CUDA_ARCHS=75;86;89' .

echo "==> tag + push ${IMAGE_URI}"
docker tag gutibm:cuda "${IMAGE_URI}"
docker push "${IMAGE_URI}"

echo "OK: pushed ${IMAGE_URI}"

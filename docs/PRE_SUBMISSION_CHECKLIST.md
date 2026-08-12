# Campaign pre-submission checklist

Run these checks before submitting a long GutIBM campaign. Local CTest is
necessary but cannot establish production MPI, GPU, AWS, or campaign-duration
behavior.

## Runtime and hardware

- [ ] Run a real graceful-stop test with `mpirun -np 2`.
- [ ] Run a real graceful-stop test with `mpirun -np 4`.
- [ ] Confirm the intended GPU is visible (`nvidia-smi` or equivalent) and run
      a physical-GPU smoke test. A GPU-labelled CTest pass in a non-CUDA build
      is not GPU coverage.
- [ ] Confirm the target AWS account, region, subnet/AZ capacity, instance
      type, quota, IAM roles, ECR image, and Batch compute environment.

## Spot interruption and restart

- [ ] Exercise the interruption path end to end:
      rank termination → closed restart → immutable upload →
      `latest.json` only after upload → status metadata.
- [ ] Confirm the closed restart is readable and contains the required agents,
      grid, summary, event, nutrient-flux, and provenance evidence.
- [ ] Confirm a resumed run starts with fresh interval counters and preserved
      cumulative counters.

## Campaign scale

- [ ] Run a campaign-duration or representative multi-hour throughput test.
- [ ] Check checkpoint cadence, upload duration, storage growth, and retained
      object sizes at the expected agent/grid scale.
- [ ] Record the binary, image, input, seed, rank count, GPU model, and AWS
      capacity details with the campaign submission.

This checklist does not replace short deterministic CI tests. CI cannot
structurally prove physical GPU execution without a visible device, arbitrary
MPI rank-count behavior, AWS Batch/S3 semantics, Spot timing, or multi-hour
throughput and storage behavior.

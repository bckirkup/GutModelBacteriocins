#include "chemical_field_gpu.h"
#include "chemical_field.h"
#include "cuda_aware_mpi.h"
#include "diffusion_gpu.h"
#include "domain.h"
#include "dispatch.h"
#include "gpu_kernels.h"
#include "device_memory.h"

#include <cstdlib>

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

void ChemicalFieldGpu::init(ChemicalField& field) {
  active_ = gpu_runtime_enabled();
  nspec_ = field.num_species();
  ncells_ = field.ncells();
  if (!active_ || nspec_ <= 0 || ncells_ <= 0) return;

  d_conc_.resize(static_cast<size_t>(nspec_));
  d_reac_.resize(static_cast<size_t>(nspec_));
  for (Int s = 0; s < nspec_; ++s) {
    d_conc_[static_cast<size_t>(s)].allocate(static_cast<size_t>(ncells_));
    d_reac_[static_cast<size_t>(s)].allocate(static_cast<size_t>(ncells_));
  }

  std::vector<double> bc(static_cast<size_t>(nspec_));
  for (Int s = 0; s < nspec_; ++s) {
    bc[static_cast<size_t>(s)] = field.spec(s).boundary_conc;
  }
  d_boundary_conc_.upload(bc);
  sync_to_device(field);
}

void ChemicalFieldGpu::sync_to_device(const ChemicalField& field) {
  sync_concentrations_to_device(field);
  sync_reactions_to_device(field);
}

void ChemicalFieldGpu::sync_to_host(ChemicalField& field) {
  sync_concentrations_to_host(field);
  sync_reactions_to_host(field);
}

void ChemicalFieldGpu::sync_concentrations_to_device(const ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    d_conc_[static_cast<size_t>(s)].upload(
        field.conc_data()[static_cast<size_t>(s)]);
  }
}

void ChemicalFieldGpu::sync_reactions_to_device(const ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<double> host(static_cast<size_t>(ncells_));
    for (Int c = 0; c < ncells_; ++c) {
      host[static_cast<size_t>(c)] = field.reac(s, c);
    }
    d_reac_[static_cast<size_t>(s)].upload(host);
  }
}

void ChemicalFieldGpu::sync_concentrations_to_host(ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<double> host(static_cast<size_t>(ncells_));
    d_conc_[static_cast<size_t>(s)].download(host);
    for (Int c = 0; c < ncells_; ++c) {
      field.conc(s, c) = host[static_cast<size_t>(c)];
    }
  }
}

void ChemicalFieldGpu::sync_reactions_to_host(ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<double> host(static_cast<size_t>(ncells_));
    d_reac_[static_cast<size_t>(s)].download(host);
    for (Int c = 0; c < ncells_; ++c) {
      field.reac(s, c) = host[static_cast<size_t>(c)];
    }
  }
}

void ChemicalFieldGpu::sync_species_concentrations_to_host(ChemicalField& field,
                                                           Int spec) {
  if (!active_ || spec < 0 || spec >= nspec_) return;
  std::vector<double> host(static_cast<size_t>(ncells_));
  d_conc_[static_cast<size_t>(spec)].download(host);
  for (Int c = 0; c < ncells_; ++c) {
    field.conc(spec, c) = host[static_cast<size_t>(c)];
  }
}

void ChemicalFieldGpu::sync_species_concentrations_to_device(
    const ChemicalField& field, Int spec) {
  if (!active_ || spec < 0 || spec >= nspec_) return;
  d_conc_[static_cast<size_t>(spec)].upload(
      field.conc_data()[static_cast<size_t>(spec)]);
}

void ChemicalFieldGpu::zero_reactions_on_device() {
#ifndef GUTIBM_CUDA
  return;
#else
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    cudaMemset(d_reac_[static_cast<size_t>(s)].data(), 0,
               static_cast<size_t>(ncells_) * sizeof(double));
  }
#endif
}

bool ChemicalFieldGpu::apply_reactions(double dt, const Domain& domain) {
  (void)domain;
#ifndef GUTIBM_CUDA
  (void)dt;
  return false;
#else
  if (!active_) return false;
  for (Int s = 0; s < nspec_; ++s) {
    gpu::launch_field_update_kernel(
        d_conc_[static_cast<size_t>(s)].data(),
        d_reac_[static_cast<size_t>(s)].data(),
        ncells_, 1, dt, gpu_compute_stream());
  }
  gpu_sync_compute();
  gpu_check_error("field_update_kernel");
  return true;
#endif
}

bool ChemicalFieldGpu::apply_diffusion(const Domain& domain,
                                       const ChemicalField& field,
                                       Real dt) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)field;
  (void)dt;
  return false;
#else
  if (!active_) return false;
  if (!gpu_diffusion_line_lengths_supported(domain)) return false;

  bool applied = false;
  for (Int s = 0; s < nspec_; ++s) {
    if (gpu_apply_species_diffusion_device(
            domain, field.spec(s), d_conc_[static_cast<size_t>(s)].data(), dt)) {
      applied = true;
    }
  }

  if (applied) {
    gpu_sync_compute();
    gpu_check_error("ChemicalFieldGpu::apply_diffusion");
  }
  return applied;
#endif
}

bool ChemicalFieldGpu::apply_boundaries(const Domain& domain,
                                        const ChemicalField& field) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)field;
  return false;
#else
  if (!active_) return false;

  const int nx = domain.nx();
  const int ny = domain.ny();
  const int nz = domain.nz();

  for (Int s = 0; s < nspec_; ++s) {
    const ChemicalSpec& spec = field.spec(s);
    double* d_conc = d_conc_[static_cast<size_t>(s)].data();
    gpu::launch_set_epithelial_boundary(
        d_conc, nx, ny, spec.boundary_conc, gpu_compute_stream());
    if (!spec.diffusion_enabled && nz >= 2) {
      gpu::launch_set_luminal_neumann(d_conc, nx, ny, nz, gpu_compute_stream());
    }
  }

  gpu_sync_compute();
  gpu_check_error("ChemicalFieldGpu::apply_boundaries");
  return true;
#endif
}

double* ChemicalFieldGpu::conc_device(Int spec) {
  if (!active_ || spec < 0 || spec >= nspec_) return nullptr;
  return d_conc_[static_cast<size_t>(spec)].data();
}

const double* ChemicalFieldGpu::conc_device(Int spec) const {
  if (!active_ || spec < 0 || spec >= nspec_) return nullptr;
  return d_conc_[static_cast<size_t>(spec)].data();
}

double* ChemicalFieldGpu::reac_device(Int spec) {
  if (!active_ || spec < 0 || spec >= nspec_) return nullptr;
  return d_reac_[static_cast<size_t>(spec)].data();
}

const double* ChemicalFieldGpu::reac_device(Int spec) const {
  if (!active_ || spec < 0 || spec >= nspec_) return nullptr;
  return d_reac_[static_cast<size_t>(spec)].data();
}

bool ChemicalFieldGpu::try_sum_reactions_on_device(ChemicalField& field) {
#ifndef GUTIBM_CUDA
  (void)field;
  return false;
#else
  if (!active_) return false;

#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return false;

  int ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  if (ranks <= 1) return false;
#endif

  const char* env = std::getenv("GUTIBM_CUDA_AWARE_MPI");
  if (env == nullptr || (env[0] != '1' && env[0] != 't' && env[0] != 'T')) {
    return false;
  }
  if (!cuda_aware_mpi_runtime_available()) {
    return false;
  }

  for (Int s = 0; s < nspec_; ++s) {
    double* d_reac = reac_device(s);
    if (d_reac == nullptr) return false;
#ifdef GUTIBM_MPI
    MPI_Allreduce(MPI_IN_PLACE, d_reac, ncells_, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
  }

  sync_reactions_to_host(field);
  return true;
#endif
}

}  // namespace gutibm

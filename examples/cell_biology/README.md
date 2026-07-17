# Cell Biology Example (Spec 3)

Short run demonstrating **Fur-regulated iron receptors**, **contact-dependent
inhibition (CDI)**, and **motility** on a compact domain.

## Enabled mechanisms

| Feature | Config keys | Fix module |
|---------|-------------|------------|
| Fur iron regulation | `fur.enabled` (default on in this example) | `fix_metabolism` |
| CDI killing | `cdi.enabled`, per-strain `cdi_type` / `cdi_immunity` | `fix_cdi` |
| Run-and-reverse motility | `motility.enabled` (+ Spec 10v2 aerotaxis/energy taxis defaults) | `fix_motility` |

Strain 1 is CDI+ with ColE1/ColB plasmids; strain 2 is CDI− immigrant.

## Run

```bash
cd build
mpirun -np 1 ./gut_ibm ../examples/cell_biology/input.json
```

See [docs/MECHANISMS.md](../../docs/MECHANISMS.md) §§7–9 and [docs/PARAMETERS.md](../../docs/PARAMETERS.md) for tunables.

# Boltz-2 run blocked on GPU allocation — P100 (sm_60)

**3 September 2026.** Four kernel versions, ~30 minutes of a sub-8-hour weekly quota. Each
failed strictly further than the last; the current blocker is hardware, not code.

| v | reached | failure |
|---|---|---|
| 1 | install | `FileNotFoundError` — bundle files nest a directory below the mount the launcher reports |
| 2 | bundle loaded (298 aa, 751 ligands) | `python -m boltz` — boltz ships a console script, not a runnable package |
| 3 | weights + CCD downloaded, 4 MSAs built, prediction started | `CUDA error: no kernel image is available for execution on the device` |
| 4 | accelerator check (< 1 min) | aborted deliberately: card is sm_60 |

## The blocker

    Tesla P100-PCIE-16GB   sm_60   torch 2.10.0+cu128
    current PyTorch supports sm_70 sm_75 sm_80 sm_86 sm_90 sm_100 sm_120

Boltz-2 installs a cu128 torch build with no kernels for Pascal. The P100 cannot run it.

## The push API cannot choose the card — now confirmed, not assumed

`acceleratorType: "nvidiaTeslaT4"` was **accepted by `/kernels/push` without error and silently
ignored**; v4 still received a P100. A rejected push costs no GPU time, so this was free to
test, and it converts a standing assumption into a measured fact.

## Two ways forward

1. **Change the notebook accelerator to T4 (sm_75) in the Kaggle UI.** A few seconds of manual
   action, no quota cost, and the timing probe then measures on the card the run will use.
2. **Downgrade torch to a cu121 build (~2.5.x), which still carries sm_60 kernels.** Doable in
   principle but it is a ~2.5 GB download inside the run, risks dependency conflicts with
   boltz's requirements and Kaggle's preinstalled stack, and any timing measured on a P100
   would not represent the hardware a real run would use. It spends scarce quota on the less
   likely path.

Option 1 is recommended. Nothing further should be pushed until the accelerator is changed:
v4's abort is now the correct behaviour and every re-push will reproduce it in under a minute.

## Not yet measured

The pre-registered probe has never run. **No Boltz-2 number exists**, so nothing about the
hypothesis is settled in either direction. The compute-abort ceiling of 6.0 GPU-hours remains
untested, and the ~20 s/ligand figure from the paper remains unverified on our hardware.

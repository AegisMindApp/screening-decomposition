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

## Resolution — no manual action needed

The first instinct was to ask for the accelerator to be switched to T4 in the Kaggle UI. Our own
records rule that out: **Kaggle has never supplied a T4, including on manual selection.** The
push API is equally powerless — `machineShape` was already known to be ignored, and this run
adds `acceleratorType` to that list.

The repo had already solved this for the scaling-law work: install a **cu121 torch that still
ships sm_60 kernels** before the package that would otherwise pull a cu128 build. `boltz 2.2.1`
declares only `torch>=2.2`, so **torch 2.4.1 satisfies it** and pip will not upgrade back to a
build the card cannot run.

v5 therefore: detects Pascal → installs `torch==2.4.1+cu121` → installs boltz → **re-verifies
with `torch.cuda.get_arch_list()`** that the installed build actually carries kernels for this
card, aborting in seconds if boltz's resolver dragged torch forward again.

## Not yet measured

The pre-registered probe has never run. **No Boltz-2 number exists**, so nothing about the
hypothesis is settled in either direction. The compute-abort ceiling of 6.0 GPU-hours remains
untested, and the ~20 s/ligand figure from the paper remains unverified on our hardware.

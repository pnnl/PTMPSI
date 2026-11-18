# Visualization Setup (ChimeraX and VMD)

Both ChimeraX and VMD are free for academic, government, and non-profit research, and neither is distributed by `pip` or Homebrew. Neither is part of PTM-Psi.

## ChimeraX

### Licensing and download
- ChimeraX requires accepting UCSF's non-commercial licence (`https://www.cgl.ucsf.edu/chimerax/docs/licensing.html`).
- Run `python .agents/skills/ptmpsi/scripts/get_chimerax.py` without flags first to print the licence notice (downloads nothing).
- Only after explicit user agreement, run `python .agents/skills/ptmpsi/scripts/get_chimerax.py --accept-license -o <dir>` (add `--platform centos-9` for Rocky Linux 9 / Deception).

### macOS (Apple Silicon, verified with ChimeraX 1.12)
- Open the downloaded `.dmg` and copy `ChimeraX-1.12.app` to `/Applications/`.
- ChimeraX opens `.gro` + `.xtc` (`open processed.gro`, then `open md_clean.xtc structureModel #1`).
- In the GUI command line, paste commands on one line separated by semicolons: `open processed.gro; open md_clean.xtc structureModel #1; show #1/A atoms; view #1`. Short capped peptides may have no cartoon representation, so use `show #1/A atoms` to display chain A.
- `--nogui` loads trajectories, selects PTM residues (`select :SNC`), and saves `.cxs` sessions, but cannot save images (`Unable to save images because OpenGL rendering is not available`). Use the GUI for image export on macOS.

### PNNL Deception (Rocky Linux 9.5, verified in job `1573564` on `dc068`)
- `chimerax` is not a cluster module. Apptainer requires `module load apptainer`; without it the executable is absent from `PATH`. Container execution in batch jobs has failed with a `squashfuse`/`fusermount` error even after loading the module. Use the unpacked RPM for the verified native offscreen workflow below. `/usr/bin/rpm2cpio` and `/usr/bin/cpio` are available.
- Unpack the `centos-9` `.rpm` into project storage without root:
  ```bash
  mkdir -p /vast/projects/<project>/<user>/chimerax-1.12-unpack
  cd /vast/projects/<project>/<user>/chimerax-1.12-unpack
  /usr/bin/rpm2cpio /path/to/ucsf-chimerax-1.12-1.el9.x86_64.rpm > pkg.cpio
  /usr/bin/cpio -idm --quiet < pkg.cpio
  rm -f pkg.cpio
  ```
- Rocky Linux 9.5's `/lib64/libcrypto.so.3` lacks `OPENSSL_3.3.0` required by ChimeraX 1.12. Preload OpenSSL from the cluster's `miniforge25.3.0` module before launching `--offscreen --nogui`:
  ```bash
  export LD_LIBRARY_PATH="/share/apps/python/miniforge25.3.0/lib:${LD_LIBRARY_PATH:-}"
  export LD_PRELOAD="/share/apps/python/miniforge25.3.0/lib/libcrypto.so.3:/share/apps/python/miniforge25.3.0/lib/libssl.so.3"
  /vast/projects/<project>/<user>/chimerax-1.12-unpack/usr/libexec/UCSF-ChimeraX/bin/ChimeraX \
    --offscreen --nogui --exit --script render.cxc
  ```

## VMD

### macOS (headless `vmd-python`, verified with `vmd-python` 3.1.7 / VMD 1.9.4a26)
- The VMD GUI installer requires a registered UIUC login download.
- For headless Python analysis/scripting without touching your shell rc files or `ptmpsi-venv`, install `micromamba` (`brew install micromamba`) and create an isolated environment:
  ```bash
  MAMBA_ROOT_PREFIX=~/.local/share/vmd-mamba micromamba create -y --no-rc -r ~/.local/share/vmd-mamba -n vmd -c conda-forge vmd-python python=3.12
  ```
- Run with `~/.local/share/vmd-mamba/envs/vmd/bin/python`. It loads `.gro` + `.xtc`, selects PTM residues (`atomsel("resname SNC")`), computes measurements, and writes external `Tachyon` scene `.dat` files (`TachyonInternal` is not compiled in).

### PNNL Deception (verified in job `1573323` on `dc068`)
- Load `module load vmd/2.0.0a7` (`vmd/1.9.4a57` is also present).
- Run `vmd -dispdev text -e script.tcl` in a batch job. `render TachyonInternal frame0.tga` renders headless on CPU nodes, and `/usr/bin/ffmpeg -y -i frame0.tga frame0.png` converts to PNG.

## Browser-compatible movies

Encode rendered frames as H.264 with an 8-bit YUV pixel format and move MP4 metadata
to the beginning of the file so that browsers can play the movie:

```bash
ffmpeg -y -framerate 20 -i frame%04d.png -c:v libx264 -pix_fmt yuv420p -movflags +faststart movie.mp4
```

### Making a trajectory movie from a Deception run (one Trp-cage WT vs SNC run)

Render on Deception, encode on the laptop:
- Deception's `/usr/bin/ffmpeg` has no `libx264` ("Unknown encoder"), so the command above fails there. Use it only for `.tga` to `.png`, then encode with a local ffmpeg.
- `dcx.py` has no download command and `read-files` is capped at 4 MiB. In a small sbatch job, `tar` the PNG frames, `base64 -w 0` the tar, `split -b 3000000` it, and record `sha256sum` of the tar. Fetch each part with `read-files`, join, decode, and check the hash locally.
- Make the molecules whole before rendering: re-extract the protein with `trjconv` only, then `MDAnalysis.lib.mdamath.make_whole`, fit on the backbone, and rotate every system with the same matrix so that they face the same way. The `gmx` nojump and fit pipeline left bad frames (see SKILL.md).
- For a side-by-side video, render each system separately, then join the frames and add labels and a timestamp locally (Pillow), then encode. Keep the header strip and the timestamp in separate areas so one cannot overwrite the other.

VMD 2.0.0a7 headless (`vmd -dispdev text -eofexit -e script.tcl < /dev/null`):
- Do not use `display rendermode GLSL` (error) or `molinfo ... get scale` (invalid attribute).
- `display resize` before a `TachyonInternal` render gave blank images. Leave the default size.
- A cartoon of a PDB with no `HELIX` records draws as a plain tube. Add `HELIX` records from DSSP.
- VMD does not treat `SNC` as a protein residue. For drawing only, write a PDB with `SNC` renamed `CYS`.
- A Tcl error ends the run silently. Wrap the body in `catch` and print the error.
- "Duplicate resname" errors at start-up are harmless.

## Windows (stub, not tested)
- ChimeraX provides a `windows` installer via `get_chimerax.py --platform windows`.
- VMD provides a Windows installer via UIUC's registered download page.
- GROMACS workflows on Windows are expected to use WSL (unverified).

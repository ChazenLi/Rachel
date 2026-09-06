# Rachel-beta Installation

Release: `Rachel-beta-1.2-20260831_120000`

This ZIP is the public runtime distribution of **Rachel-beta**. It contains the
Codex Skill and the built-in base knowledge required for normal operation. It
does not contain the internal development workspace, tests, research records,
system-overview packet, private knowledge packs, or staging data.

## 1. Prepare the environment

```powershell
conda env create -f environment.yml
conda activate rachel-v2
```

The `rachel-v2` Conda environment name is a compatibility identifier; the
product name is Rachel-beta.

The public package is verified on Windows with Python 3.11 and RDKit. Other
operating systems are not claimed as verified by this distribution.

## 2. Verify the extracted release

From the extracted release root:

```powershell
python VERIFY_RELEASE.py
```

Continue only when the JSON output contains `"ok": true`. The verifier checks
the package checksums, runtime imports, base profile, reaction-site disclosure,
validation closure, and a small non-mutating route smoke.

## 3. Install the Skill

Use `$env:CODEX_HOME` when it is configured. Otherwise Codex defaults to
`$HOME\.codex` on Windows.

```powershell
$CodexHome = if ($env:CODEX_HOME) { $env:CODEX_HOME } else { Join-Path $HOME '.codex' }
$SkillTarget = Join-Path $CodexHome 'skills\Rachel'
New-Item -ItemType Directory -Force -Path (Split-Path $SkillTarget) | Out-Null
Copy-Item -Recurse -Force '.\Rachel' $SkillTarget
```

Restart Codex or open a new task so the Skill catalog is refreshed.

## 4. License

Rachel-beta is distributed under CC BY-NC-ND 4.0 International. Attribution,
non-commercial use, and no-derivatives conditions apply. See `LICENSE`.

# Rachel-beta Release Notes

Release: `Rachel-beta-1.2-20260831_120000`
Date: 2026-08-31T15:32:59.793196+08:00

## Public distribution boundary

This release is a deliberately narrow runtime ZIP. It includes the executable
Rachel Skill, built-in chemistry and validation resources, base knowledge, a
minimal public README, installation instructions, checksums, and the packaged
verifier. Internal tests, research and design documents, walkthrough sessions,
caches, staging workspaces, private team/project knowledge packs, and local
development paths are excluded.

## Included beta capabilities

- stateful route construction and explicit `try_action -> validation -> commit`;
- site-first action discovery and complete custom-action registration;
- canonical `rachel.validation.v2` findings and route-digest-bound audit;
- opt-in current-node Scouting through `scout_node` and `scout_record`;
- immutable built-in base profile with optional `base -> team -> project`
  profile composition supplied by the user;
- controlled outcomes and inactive knowledge drafts in the runtime contract;
- report, export, visualization, and clean-runtime verification paths.

## Verification boundary

The package is verified on Windows with Python 3.11 and RDKit. The maintainer
test suite and internal release checks are run outside this ZIP. This release
does not claim that arbitrary chemical routes are experimentally validated,
that `clear` means experimental success, or that external private-pack
publication has a complete independent end-to-end regression module.

## License

CC BY-NC-ND 4.0 International applies to this distribution. See `LICENSE` for
the attribution, non-commercial, and no-derivatives conditions.

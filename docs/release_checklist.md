# Release and archiving checklist

Software versions and DOIs are separate from the dataset's: the processed
temperature dataset is version 1.0.0 with DOI 10.5281/zenodo.21717874
(concept DOI 10.5281/zenodo.21717873). The software release described here
is `v1.0.0` of this repository and gets its own DOI when archived.

## Steps for a release

1. Update `CHANGELOG.md`, the version in `pyproject.toml`,
   `src/scorch/__init__.py`, `CITATION.cff`, `codemeta.json` and
   `.zenodo.json` (all `1.0.0` for this release) and `date-released`.
2. Run the tests, including the study regression test with the archived
   input (`SCORCH_STUDY_INPUT=input/tmax_processed.nc python -m pytest`), and
   record the result in `docs/reproduction_record.md`.
3. Commit and push to `main`; wait for the GitHub Actions workflow to pass.
4. Create the GitHub release (tag `v1.0.0` on the tested commit) with the
   changelog entry as release notes, for example
   `gh release create v1.0.0 --target <commit> --title "SCORCH v1.0.0" --notes-file <notes>`.
5. Archive on Zenodo. Either enable the repository in the Zenodo GitHub
   integration (https://zenodo.org/account/settings/github/) **before**
   creating the release, so that Zenodo archives the release automatically
   and reads the metadata from `.zenodo.json`; or upload the release archive
   (`Source code (zip)` of the release) to a new Zenodo software record by
   hand, using the metadata of `.zenodo.json`. Zenodo mints the DOI on
   publication; it cannot be undone, so check the metadata first.
6. Put the software DOI into `CITATION.cff` (`identifiers:` with
   `type: doi`), into the README badge and into `CHANGELOG.md`, and add the
   relation "is supplement to" / "references" between the software record
   and the dataset record on Zenodo. Commit that small follow-up.
7. For the article, fill the code metadata table from
   `softwarex_checklist.md` and cite the software DOI and the dataset DOI.

The Zenodo steps need the account owner's login; they were not done for
v1.0.0 at the time of the release commit (see `reproduction_record.md`).

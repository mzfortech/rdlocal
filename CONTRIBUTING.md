# Contributing to rdlocal

Thanks for your interest in improving `rdlocal`. Bug reports, documentation fixes, and carefully validated statistical contributions are welcome.

## Development setup

Use Python 3.11 or newer. From the repository root, create and activate a virtual environment, then install the package and development dependencies:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Run the regular test suite with:

```sh
pytest
```

The Monte Carlo coverage tests are marked `slow`; run them separately with `pytest -m slow` when changing statistical calculations.

To build the documentation locally, install `docs/requirements.txt` and run `python -m sphinx -W --keep-going -b html docs docs/_build/html`.

## Clean-room research rules

`rdlocal` is an independent implementation based on published academic papers. Keep contributions within that boundary:

- Use papers and their supplements as the specification. Cite the relevant paper, section, equation, or page in the code or pull request when changing a statistical method.
- Do not consult or copy source code, docstrings, help files, manuals, vignettes, READMEs, or other repository files from `rdpackages` software. The full project rules and source record are in [`CLEANROOM.md`](CLEANROOM.md).
- Existing reference numbers were produced by black-box comparisons. New comparisons may use other implementations only as black boxes; record the inputs, call arguments, and returned numeric results, never their code or explanatory text.
- Use synthetic or otherwise redistributable data in tests. Do not commit the Senate vote-margin data or other restricted data.
- When a paper leaves a step unspecified, describe the chosen interpretation and its limitations instead of presenting it as settled. Update `CLEANROOM.md` when a change adds or resolves a methodological deviation.

## Tests and pull requests

Submit changes through a pull request; do not push directly to `main`. A pull request needs at least one approval from someone other than its author, all required CI checks must pass, and open review conversations must be resolved before merge. If you update a pull request after approval, request another review.

Add a focused test for bug fixes and statistical changes. Prefer deterministic inputs and small examples; keep long simulation checks under the `slow` marker. Run the regular suite before opening a pull request, and run the slow suite when changing coverage or simulation behavior.

In the pull request description, explain the problem, summarize the change, list the tests run, and cite research supporting any statistical change. Keep changes focused and update the README or methodology notes when user-visible behavior or limitations change.

Contributions are made under the project's BSD-3-Clause license, as shown in `LICENSE`.

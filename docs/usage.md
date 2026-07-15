# Usage

Use a CSV with unique `sample_id`, `idat_red`, and `idat_green` columns. Existing red/green files are required, may be gzipped, and must form a matching basename pair. Optional metadata are retained verbatim and can be used in `--design`.

The default SeSAMe preparation code is `QCDPB` and is passed explicitly to `prepSesame()`. The container records SeSAMe and R versions per sample. Provide a custom `--sesame_prep_code` only after validating it against the container SeSAMe release.

Cell composition is disabled by default; no bundled reference is represented as EPICv2-validated. Include known batch factors in the DMP formula rather than destructively correcting the principal matrices.

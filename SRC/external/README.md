# External vendored libraries

## toml-f
Fortran TOML parser. License: Apache-2.0 OR MIT.
Source: https://github.com/toml-f/toml-f, v0.5.0.
See `toml-f/VENDORED_FROM.txt` for details and local patches.

ecalj uses toml-f to read `ctrlg.<sname>.toml` (m_ctrl_toml_loader for the lmf side, m_GWinput for the GW side).

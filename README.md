# ecalj

First-principles electronic-structure package (PMT = APW + MTO basis; LDA/GGA and QSGW).

This package is meant to be read by Claude (Claude Code), which then explains it to people: open the repository in Claude Code
and ask. Claude starts from `CLAUDE.md` (it reads `MD/ecaljclaude.md`); the notes are in `MD/` (`MD/handover.md` for a session without memory), the changes in `Changes.txt`.
The manual for people is [ecaljdoc](https://ecalj.github.io/ecaljdoc/); its source is the directory `ecaljdoc/` of this repository
(`ecaljdoc/manual/`, `ecaljdoc/theory/`; the development notes `MD/` are not on the site). Each sample in `Samples/` has a README
that points to its section of the manual (`ecaljdoc/manual/samples.md` lists them). `python3 TOOLS/doclinks.py` checks the links.
